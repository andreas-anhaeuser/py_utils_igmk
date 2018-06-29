#!/usr/bin/python
"""Performance info for loops.

    See docstring of Chronometer for more info.

    Author
    ------
    Andreas Anhaeuser (AA) <andreas.anhaeuser@posteo.net>
    Institute for Geophysics and Meteorology
    University of Cologne, Germany

    History
    -------
    2015       (AA): Created
    2017-01-25 (AA): Changed layout
    2018-01-25 (AA): Updated docstring.
"""
import os
import datetime as dt
import warnings

import numpy as np

import string_utils

###################################################
# CONSTANTS                                       #
###################################################
# column widths
_col_width = [11] + [20] * 3

# time formats
_tfmt = '%Y-%m-%d %H:%M:%S'
_tfmt_days = '%d %H:%M:%S'
_tfmt_hours = '%d %H:%M:%S'
_tfmt_minutes = '%d %H:%M:%S'

# foreground colors
_ENDC = string_utils._ENDC
_BOLD = string_utils._BOLD
_BLUE = string_utils._BLUE
_YELLOW = string_utils._YELLOW
_RED = string_utils._RED
_GREEN = string_utils._GREEN

# motions
_UPWARD = '\033[1A'
_DOWNWARD = '\033[1B'
_FORWARD = '\033[1C'
_BACKWARD = '\033[1D'

_CLRSCR = '\033[2J'     # clear screen, move to (0, 0)
_CLEARLINE  = '\033[K'  # erase to end of line
_SAVECRS = '\033[s'     # save curser position
_RESTORECRS = '\033[u'  # restore curser position

###################################################
# CLASS                                           #
###################################################
class Chronometer(object):
    """A performance info screen for expensive loops.
    
        This is used to show performance info on screen during execution of
        time-consuming loops. The (approximate) number of iterations must be
        known beforehand.

        Howto
        -----
        1. Create an instance of this class before (outside) your loop.
        2. At each iteration call the loop() method.
        3. Call show() at the point where you want the status message to be
           updated.
           - often, loop() and show() will be called on the same place of your
             code. In this case, you can call loop_and_show() which has the
             same effect.
        4. If you want to print other messages on screen, use the method
           issue(). Otherwise (i. e. if you use `print` instead) the screen
           text will get distorted.
        5. After the loop it is recommended to call resumee(). This will update
           the status screen is
           updated a last time.

        Caution
        -------
        If you use the built-in `print` to issue text while Chronometer is
        active, the screen text will get distorted. Use Chronometer.issue()
        instead.

        Example code
        ------------
        >>> N = 10**6  # number of iterations
        >>> chrono = Chronometer(N, header='my meaningful header')
        >>> for n in range(N):
        ...    chrono.loop()
        ...    do_something_easy()
        ...    
        ...    # special case: jump to next loop before doing anything
        ...    # expensive
        ...    if something_extraordinary():
        ...        chrono.decrease_total_count()
        ...        continue
        ...
        ...    # the real expensive part
        ...    do_something_tough()
        ...
        ...    # update status screen
        ...    chrono.show()
        ...
        ...    # use `issue` instead of `print` while Chronometer is active
        ...    chrono.issue(some_message_string)
        >>>
        >>> # Update for the last time (and slightly modify appearance)
        >>> chrono.resumee()
    """
    ###################################################
    # INIT                                            #
    ###################################################
    def __init__(self, total_count, time_step=0.02, header='',
            show_message_times=True,):
        """Initialize.

            Parameters
            ----------
            total_count : int
                estimated number of iterations
            time_step: float, optional
                (seconds) time interval in which messages should be printed on
                screen. Default: 0.02 (fast enough for the human eye, slow
                enough to not bother the machine very much).
            header : str, optional
                header for the message that is printed on screen.
            show_message_times : bool, optional
                Whether to prefix issue() messages by time stamp.
        """
        now = dt.datetime.now()

        self.total_count = total_count
        self.count = 0

        self.header = header

        self.time_step = time_step
        self.start = now

        self.last_active = now
        self.Nlines_last_message = 0
        self.time_of_last_message = now
        self.show_message_times = show_message_times

        warnings.showwarning = self.showwarning

    ###################################################
    # HELPER FUNCTIONS                                #
    ###################################################
    def speed_string(self):
        speed = self.speed()
        if speed == 0:
            return '---'
        x = speed
        units = 's'
        # minutes
        if x < 1:
            x *= 60
            units = 'min'
        # hours
        if x < 1:
            x *= 60
            units = 'h'
        # days
        if x < 1:
            x *= 24
            units = 'd'
        return '%s loops/%s' % (string_utils.human_format(x), units)

    def inverse_speed_string(self):
        speed = self.speed()
        if speed == 0:
            return '---'
        inverse_speed = 1 / speed
        return short_time_string(inverse_speed, sep=' ') + '/loop'

    def now_string(self):
        now = dt.datetime.now()
        return now.strftime(_tfmt)

    def start_string(self):
        start = self.start
        return start.strftime(_tfmt)

    def end_string(self):
        now = dt.datetime.now()
        time_todo = self.time_todo()
        if time_todo is None:
            return '---'
        end = now + time_todo
        return end.strftime(_tfmt)

    def prefix_for_issue(self):
        if not self.show_message_times:
            return ''

        now = dt.datetime.now()
        last = self.time_of_last_message
        self.time_of_last_message = now
        seconds_since_last = (now - last).total_seconds()

        now_str = now.strftime('%H:%M:%S')
        # time_since_last_str = '{:3.0f}s'.format(seconds_since_last)
        time_since_last_str = short_time_string_no_frac(
                seconds_since_last, 2).rjust(3)
        prefix = '%s [%s] ' % (now_str, time_since_last_str)

        return _BLUE + prefix + _ENDC

    def build_line(self, words, colors=[_BOLD, None, None, None]):
        line = ''
        for i, word in enumerate(words):
            width = _col_width[i] - 1
            color = colors[i]
            word = word.rjust(width) 
            if color is not None:
                word = color + word + _ENDC
            line = line + word + ' '
        line = line + '\n'
        return line

    def get_status_text(self, usermessage=None, mode='run'):
        """Return a formatted string showing the progress.

            Parameters
            ----------
            usermessage : str, optional
                Text to be shown below the progress overview. May be several
                lines.
            mode : {'run', 'resumee'}, optional
                Has impact on how parts of the text are color highlighted.
                Default: 'run'

            Returns
            -------

            Author
            ------
            Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
            Institute for Geophysics and Meteorology
            University of Cologne, Germany

            History
            -------
            2015       (AA): Created
            2017-04-06 (AA): Added parameter `mode`
        """
        ###################################################
        # INPUT CHECK                                     #
        ###################################################
        if mode is None:
            mode = 'run'
        valid_modes = ['run', 'resumee']
        assert mode in valid_modes

        ###################################################
        # HEADER                                          #
        ###################################################
        # initialize
        text = ''
        words = [''] * len(_col_width) 

        # empty line
        text = text + '\n'

        # header line
        if self.header != '':
            length = sum(_col_width)
            line = _BLUE + self.header.center(length) + _ENDC + '\n\n'
            text = text + line

        ###################################################
        # SPEED AND TIME                                  #
        ###################################################
        _colors = (_BOLD, None, _BOLD, None)

        # first line
        words[0] = 'Start:'
        words[1] = self.start_string()
        words[2] = ''
        words[3] = ''
        line = self.build_line(words, colors=_colors)
        text = text + line

        # second line
        if mode == 'run':
            colors = (_BOLD,  _GREEN, _BOLD, None)
        elif mode == 'resumee':
            colors = (_BOLD,  None, _BOLD, None)
        words[0] = 'Now:'
        words[1] = self.now_string()
        words[2] = 'Speed:'
        words[3] = self.speed_string()
        line = self.build_line(words, colors=colors)
        text = text + line

        # third line
        words[0] = 'End:'
        words[1] = self.end_string()
        words[2] = 'Inverse speed:'
        words[3] = self.inverse_speed_string()
        if mode == 'run':
            colors = (_BOLD,  _YELLOW, _BOLD, None)
        elif mode == 'resumee':
            colors = (_BOLD,  _GREEN, _BOLD, None)
        line = self.build_line(words, colors=colors)
        text = text + line

        # empty line
        line = '\n'
        text = text + line

        # column headers
        words = ['', 'Total', 'Done', 'Remaining']
        colors = (_BOLD,) * len(words)
        line = self.build_line(words, colors=colors)
        text = text + _BOLD + line + _ENDC

        _colors = (_BOLD, None, None, None)

        # time
        words[0] = 'Time:'
        words[1] = time_string(self.time_total())
        words[2] = time_string(self.time_done())
        words[3] = time_string(self.time_todo())
        if mode == 'run':
            colors = (_BOLD,  None, None, _YELLOW)
        elif mode == 'resumee':
            colors = (_BOLD,  _GREEN, None, None)
        line = self.build_line(words, colors=colors)
        text = text + line

        # count line
        words[0] = 'Loops:'
        words[1] = count_string(self.total_count)
        words[2] = count_string(np.floor(self.count))
        words[3] = count_string(np.ceil(self.total_count - self.count))
        line = self.build_line(words, colors=_colors)
        text = text + line

        # fraction line
        words[0] = 'Fraction:'
        words[1] = '100 %'
        words[2] = string_utils.percentage_string(self.fraction_done())
        words[3] = string_utils.percentage_string(self.fraction_todo())
        line = self.build_line(words, colors=_colors)
        text = text + line

        # bar
        bar_width = sum(_col_width[1:])
        fraction_done = self.fraction_done()
        bar = string_utils.progress_bar(fraction_done, bar_width)
        text = text + ' ' * _col_width[0] + bar + '\n'

        # usermessage
        if usermessage is not None:
            line = '\n' + usermessage + '\n'
            text = text + line

        return text

    def update_screen(self, text):
        Nlines = self.Nlines_last_message
        for nline in range(Nlines):
            print(_UPWARD + _CLEARLINE + _UPWARD)
        print(text)

        Nlines = text.count('\n')
        self.Nlines_last_message = Nlines + 1

    def fraction_done(self):
        total_count = self.total_count
        if total_count > 0:
            return 1. * self.count / self.total_count
        else:
            return 1.

    def fraction_todo(self):
        return 1. - self.fraction_done()

    def time_done(self):
        """Return a dt.timedelta."""
        return dt.datetime.now() - self.start

    def time_todo(self):
        """Return a dt.timedelta."""
        if self.count == self.total_count:
            return dt.timedelta()

        time_total = self.time_total()
        if time_total is None:
            return None
        return time_total - self.time_done()

    def time_total(self):
        """Return a dt.timedelta."""
        seconds_done = self.time_done().total_seconds()
        fraction_done = self.fraction_done()
        if fraction_done == 0:
            return None
        seconds_total = seconds_done / fraction_done
        return dt.timedelta(seconds=seconds_total)

    def speed(self):
        seconds_done = self.time_done().total_seconds()
        return self.count / seconds_done
        
    ###################################################
    # USER FUNCTIONS                                  #
    ###################################################
    def loop(self):
        """Equivalent to self.increase_count(1)."""
        self.count += 1

    def show(self, usermessage=None, force=None, mode=None):
        if force is None:
            if usermessage is None:
                force = False
            else:
                force = True
        
        now = dt.datetime.now()
        time_inactive = now - self.last_active
        if not force:
            if time_inactive.total_seconds() < self.time_step:
                return
        text = self.get_status_text(usermessage=usermessage, mode=mode)
        self.update_screen(text)
        self.last_active = now

    def loop_and_show(self, usermessage=None, force=None):
        self.loop()
        self.show(usermessage=usermessage, force=force)

    def issue(self, text):
        prefix = self.prefix_for_issue()
        self.update_screen(prefix + text)
        self.Nlines_last_message = 0
        self.show(force=True)

    def resumee(self, usermessage=None):
        self.show(force=True, usermessage=usermessage,
                mode='resumee',)

    def warning(self, message, prefix='WARNING: '):
        """Issue a warning."""
        self.issue(_YELLOW + prefix+ _ENDC + message )

    def showwarning(self, *args, **kwargs):
        """Overwrite default implementation."""
        prefix = ' ' * 24
        warning = args[0]
        type = args[1]
        where = args[2]
        lineno = args[3]
        message = warning.message
        text = (message + ' in\n' +
                prefix + str(lineno) + ':' + where + '\n' +
                prefix + '(Type: ' + str(type) + ')')
        self.warning(text)

    def debug_warning(self, module_name):
        """Issue a DEBUG-mode warning."""
        text = (_BOLD + _RED + 'DEBUG' + _ENDC +
                 '-mode in ' +
                 _BOLD + module_name + _ENDC)
        self.warning(text)

    ###################################################
    # COUNT                                           #
    ###################################################
    def set_count(self, c):
        self.count = c

    def get_count(self):
        return self.count

    def decrease_count(self, increment=1):
        self.count -= increment

    def increase_count(self, increment=1):
        self.count += increment

    ###################################################
    # TOTAL COUNT                                     #
    ###################################################
    def set_total_count(self, c):
        self.total_count = c

    def get_total_count(self):
        return self.total_count

    def decrease_total_count(self, increment=1):
        self.total_count -= increment

    def increase_total_count(self, increment=1):
        self.total_count += increment

class PerformanceInfo(Chronometer):
    """Alias to Chronometer."""
    pass

###################################################
# HELPER FUNCTIONS                                #
###################################################
def sound(length, freq=1000):
    """Play a sine sound.

        Parameters
        ----------
        length : float
            duration in seconds
        freq : float
            frequency in Hertz

        Returns
        -------
        None

        Author
        ------
        Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
        Institute for Geophysics and Meteorology
        University of Cologne, Germany

        History
        -------
                2016-09-22 (AA): Created
    """
    os.system('play --no-show-progress --null --channels 1 synth %s sine %f' %\
            (length, freq))

def time_string(x):
    """Convert timedelta to str.

        Parameters
        ----------
        x : dt.timedelta or None

        Returns
        -------
        str

        Author
        ------
        Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
        Institute for Geophysics and Meteorology
        University of Cologne, Germany

        History
        -------
        2017-04-06 (AA): Created
    """
    ###################################################
    # SPECIAL CASES                                   #
    ###################################################
    # None
    if x is None:
        return 'unknown'

    # negative
    if x < dt.timedelta():
        return '-' + time_string(-x)
    
    ###################################################
    # REGULAR CASE                                    #
    ###################################################
    secs = int(x.total_seconds())
    return str(dt.timedelta(seconds=secs))

def short_time_string(seconds, sep=''):
    """Convert to minutes, hours, days if neccessary.

        Parameters
        ----------
        seconds : float
        sep : str, optional
            separator between number and unit. Default: ''

        Returns
        -------
        str
            something like '3.0s', '51min', etc

        Author
        ------
        Andreas Anhaeuser (AA) <anhaeus@meteo.uni-koeln.de>
        Institute for Geophysics and Meteorology
        University of Cologne, Germany

        History
        -------
        2017-03-08 (AA): Created
    """
    x = seconds
    units = 's'
    # minutes
    if abs(x) > 60:
        x /= 60
        units = 'min'
    # hours
    if abs(x) > 60:
        x /= 60
        units = 'h'
        # days
        if abs(x) > 24:
            x /= 24
            units = 'd'
    return '%s%s' % (string_utils.human_format(x, sep=sep), units)

def short_time_string_no_frac(seconds, digits=2, sep=''):
    x = int(np.round(seconds))
    units = 's'
    # minutes
    if len(str(x)) > digits:
        x /= 60
        units = 'm'
    # hours
    if len(str(x)) > digits:
        x /= 60
        units = 'h'
    # days
    if len(str(x)) > digits:
        x /= 24
        units = 'd'
    # years
    if len(str(x)) > digits:
        x /= 365
        units = 'y'
    return '%s%s%s' % (str(x), sep, units)

def count_string(count):
    if count < 0:
        return '-' + count_string(-count)

    plain = '%1.0f' % count
    len_plain = len(plain)
    nice = plain[-3:]
    pos = -3
    while pos > - len_plain:
        nice = plain[pos-3:pos] + ',' + nice
        pos -= 3
    return nice


###################################################
# TESTING                                         #
###################################################
if 1 and __name__ == '__main__':
    N = 2**16
    K = 10

    for case in (0, 1):
        if case == 0:
            Nchrono = K * N
            header = 'Chronometer'
        elif case == 1:
            Nchrono = K * N * 8/10
            header = 'Chronometer (with wrong estimate of total loops)'

        chrono = Chronometer(Nchrono, header=header)
        for k in range(K):
            chrono.issue('loop %i' % k)
            for n in range(N):
                chrono.loop()
                chrono.show()

                np.cos(np.cos(np.cos(1)))
                
        chrono.resumee('Done')
