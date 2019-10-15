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
from __future__ import print_function

# standard modules
import os
import sys
from collections import Iterable
from copy import deepcopy as copy
import datetime as dt
import inspect
import textwrap
import warnings
if sys.version_info.major < 3:
    import __builtin__ as builtin
else:
    import builtins as builtin

# PyPI modules
import numpy as np

# local modules
if __name__ == '__main__':
    import string_utils
    from timer import Timer
    import chronometer_utils as utils
else:
    from . import string_utils
    from .timer import Timer
    from . import chronometer_utils as utils

_builtin_print = copy(builtin.print)
_builtin_warning = warnings.showwarning

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
    def __init__(
            self, total_count=1, time_step=None, header='', info='',
            show_message_times=True, file=None, print_colors=None,
            item_name='loop', item_plural=None, report_level='normal',
            ):
        """Initialize.

            Parameters
            ----------
            total_count : int or Iterable
                estimated number of iterations
            time_step: float, optional
                (seconds) time interval in which messages should be printed on
                screen. Default: 0.02 (fast enough for the human eye, slow
                enough to not bother the machine very much).
            header : str, optional
                header for the message that is printed on screen.
            show_message_times : bool, optional
                Whether to prefix issue() messages by time stamp.
            item_name : str, optional
                (default: 'loop') What is counted.
            item_plural : str, optional
                (default: `unit` + '(e)s') Plural if not standard (item_name +
                (e)s).


            TODO
            ====
            Report levels
            -------------
            1 : critical
            2 : error
            3 : warning
            4 : normal
            5 : debug
        """
        self.set_output(file, print_colors, time_step)
        self.set_total_count(total_count)
        self.set_count(0)
        self.set_header(header)
        self.info = str(info)
        self.set_item_name(item_name, item_plural)

        # timers
        self.global_timer = Timer().start()
        self.last_active_timer = Timer().start()
        self.last_message_timer = Timer().start()

        self.Nlines_last_message = 0
        self.show_message_times = show_message_times

        # override builtin `print` with custom print function
        self.use_object_print()

        # override default warning handling
        warnings.showwarning = self.custom_warning

    ############################################################
    # CONTEXT MANAGEMENT                                       #
    ############################################################
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Clean up."""
        self.exit()
        self.resumee()

    def __del__(self):
        self.exit()

    ###################################################
    # USER FUNCTIONS                                  #
    ###################################################
    @classmethod
    def exit(self, usermessage=None):
        """Re-assign builtin functions for `print` and warnings."""
        builtin.print = _builtin_print
        warnings.showwarning = _builtin_warning
        return self

    def loop(self, loops=1):
        """Equivalent to self.increase_count(1)."""
        self.count += loops

    def show(self, usermessage=None, force=None, mode=None, wrap=True):
        """Update screen."""
        # set force to True of False
        if force is None:
            if self.count < 2:
                force = True
            elif usermessage is not None:
                force = True
            else:
                force = False

        # check if need to show
        now = dt.datetime.now()
        inactive_sec = self.last_active_timer.get('s')
        if inactive_sec < self.time_step and not force:
            return self

        # python2-compatibility
        python_version = sys.version_info[0]    # (int)
        if python_version < 3:
            usermessage = unicode(usermessage)

        if wrap is True:
            wrap = self.get_wrap_length()

        if wrap and (usermessage is not None):
            lines = textwrap.wrap(
                    usermessage, int(wrap), break_on_hyphens=False,
                    )
            usermessage = ''.join([line + '\n' for line in lines])

        self.last_active_timer.reset()

        text = self.get_status_text(usermessage=usermessage, mode=mode)
        self.update_screen(text)

        return self

    def loop_and_show(self, usermessage=None, force=None):
        self.loop()
        self.show(usermessage=usermessage, force=force)

    def issue(self, *args, **kwargs):
        return self.print(*args, **kwargs)

    def use_builtin_print(self):
        builtin.print = _builtin_print
        return self

    def use_object_print(self):
        builtin.print = self.print
        return self

    def print(self, text, *args, **kwargs):
        """Print time stamp and message."""
        if 'file' in kwargs:
            _builtin_print(text, *args, **kwargs)

        warnings.showwarning = _builtin_warning

        if 'sep' in kwargs:
            sep = kwargs['sep']
        else:
            sep = ' '

        if 'wrap' in kwargs:
            wrap = kwargs['wrap']
        else:
            wrap = None

        if wrap is None:
            if self.file is None:
                wrap = True
            else:
                wrap = False

        # join arguments to one string
        text = sep.join([str(arg) for arg in (text,) + args])

        # wrap text ----------------------------------
        if wrap:
            if isinstance(wrap, bool):
                wrap = self.get_wrap_length()

            lines = textwrap.wrap(
                text, wrap, break_on_hyphens=False
                )
        else:
            lines = [text]

        if len(lines) == 0:
            lines = ['']
        # --------------------------------------------

        # print --------------------------------------
        for nline, line in enumerate(lines):
            if nline == 0:
                prefix = self.prefix_for_issue()
            else:
                prefix = ' ' * self.prefix_length()
            self.update_screen(prefix + line)
            self.Nlines_last_message = 0
            self.show(force=True)
        # --------------------------------------------
            
        warnings.showwarning = self.custom_warning
        return self

    def get_wrap_length(self):
        with os.popen('stty size', 'r') as fid:
            screen_size = fid.read().split()
        wrap_length = int(screen_size[1]) - (self.prefix_length() + 1)
        return wrap_length

    def resumee(self, usermessage=None):
        """Show an overview and clean up."""
        self.show(force=True, usermessage=usermessage, mode='resumee')
        self.exit()
        return self

    def warning(self, message, prefix='WARNING: '):
        """Issue a warning."""
        if self.print_colors:
            self.print(_YELLOW + prefix + _ENDC + str(message))
        else:
            self.print(prefix + str(message))

    def custom_warning(self, *args, **kwargs):
        """Overwrite default implementation."""
        prefix = ' ' * 24
        warning = args[0]
        type = args[1]
        where = args[2]
        lineno = args[3]

        # ========== retrieve message =================== #
        # how to do this depends on major python version
        python_version = sys.version_info[0]    # (int)
        if python_version <= 2:
            message = warning.message
        else:
            message = warning.args[0]
        # ================================================ #

        text = (message + ' in\n' +
                prefix + str(lineno) + ':' + where + '\n' +
                prefix + '(Type: ' + str(type) + ')')
        self.warning(text)

    def debug_warning(self, show_warning=True, module_name=None):
        """Issue a DEBUG-mode warning.

            Parameters
            ----------
            show_warning : bool, optional
                (default: True) If False, do nothing
            module_name : str or None, optional
                if None, module name will be determined automatically
        """
        if not show_warning:
            return self

        # retrieve name of calling function
        if module_name is None:
            module_name = inspect.stack()[1][1]

        # compose message (with fancy colors and so)
        if self.print_colors:
            text = (_BOLD + _RED + 'DEBUG' + _ENDC
                     + '-mode in '
                     + _BOLD + module_name + _ENDC)
        else:
            text = 'DEBUG-mode in ' + module_name

        # issue message
        self.warning(text)

        return self

    ###################################################
    # GETTERS                                         #
    ###################################################
    def get_count(self):
        return self.count

    def get_info(self):
        return self.info

    def get_total_count(self):
        return self.total_count

    ###################################################
    # SETTERS                                         #
    ###################################################
    def set_count(self, c):
        self.count = c
        return self

    def set_file(self, *args, **kwargs):
        """Alias to set_output."""
        return self.set_output(*args, **kwargs)

    def set_output(self, filename=None, print_colors=None, time_step=None):
        """Set output file, colors, time step."""
        if print_colors is None:
            if filename is None:
                print_colors = True
            else:
                print_colors = False

        if time_step is None:
            if filename is None:
                time_step = 0.03
            else:
                time_step = 1.

        self.print_colors = print_colors
        self.time_step = time_step
        self.file = filename
        if self.file is None:
            self.file_initialized = True
        else:
            self.file_initialized = os.path.isfile(self.file)
        return self

    def set_info(self, info=''):
        """Change info string."""
        if not isinstance(info, str):
            raise TypeError('info must be str.')
        self.info = str(info)

    def set_header(self, header):
        self.header = str(header)
        return self

    def set_item_name(self, item_name=None, plural_name=None):
        if item_name is None:
            item_name = 'loop'

        if plural_name is None:
            plural_name = utils.plural(item_name)

        self.item_name = item_name
        self.item_plural = plural_name
        return self

    def set_total_count(self, c):
        if isinstance(c, Iterable):
            self.total_count = len(c)
        else:
            self.total_count = int(c)
        return self

    ###################################################
    # COUNT                                           #
    ###################################################
    def decrease_count(self, increment=1):
        self.count -= increment
        return self

    def increase_count(self, increment=1):
        self.count += increment
        return self

    ###################################################
    # TOTAL COUNT                                     #
    ###################################################
    def decrease_total_count(self, increment=1):
        self.total_count -= increment
        return self

    def increase_total_count(self, increment=1):
        self.total_count += increment
        return self

    ###################################################
    # HELPER FUNCTIONS                                #
    ###################################################
    def prefix_for_issue(self):
        """Return a colored time stamp as str."""
        if not self.show_message_times:
            return ''

        # absolute time
        now_str = dt.datetime.now().strftime('%H:%M:%S')

        # difference to last message
        diff = self.last_message_timer.get('s')
        diff_str_inner = utils.nice_time_string(diff)
        diff_str_full = '+%s' % diff_str_inner.rjust(4)
        self.last_message_timer.reset()

        # combine
        prefix = '%s [%s] ' % (diff_str_full, now_str)

        if self.print_colors:
            return _BLUE + prefix + _ENDC
        else:
            return prefix

    def prefix_length(self):
        pattern = '+0.0s [00:00:00] '
        return len(pattern)

    def build_line(self, words, colors=[_BOLD, None, None, None]):
        line = ''
        for i, word in enumerate(words):
            width = _col_width[i] - 1
            color = colors[i]
            word = word.rjust(width)
            if color is not None and self.print_colors:
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
        valid_modes = ('run', 'resumee')
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
        indent = ' ' * _col_width[0]
        header = self.header
        if mode == 'resumee':
            header = '[FINISHED] ' + header
        if header != '':
            if self.print_colors:
                line = indent + _BLUE + self.header + _ENDC + '\n'
            else:
                line = indent + self.header + '\n'
            text = text + line

        if self.info != '':
            lines = self.info.split('\n')
            for line in lines:
                if self.print_colors:
                    line = indent + line + _ENDC + '\n'
                else:
                    line = indent + line + '\n'
                text = text + line

        text = text + '\n'

        ###################################################
        # SPEED AND TIME                                  #
        ###################################################
        _colors = (_BOLD, None, _BOLD, None)

        # first line
        words[0] = 'Start:'
        words[1] = self.global_timer.get_start().strftime(_tfmt)

        words[2] = ''
        words[3] = ''
        line = self.build_line(words, colors=_colors)
        text = text + line

        # second line
        if mode == 'run':
            colors = (_BOLD, _GREEN, _BOLD, None)
        elif mode == 'resumee':
            colors = (_BOLD, None, _BOLD, None)
        words[0] = 'Now:'
        words[1] = dt.datetime.now().strftime(_tfmt)
        words[2] = 'Speed:'
        words[3] = utils.speed_string(self)
        line = self.build_line(words, colors=colors)
        text = text + line

        # third line
        words[0] = 'End:'
        words[1] = utils.end_string(self, _tfmt)
        words[2] = 'Expenditure:'
        words[3] = utils.inverse_speed_string(self)
        if mode == 'run':
            colors = (_BOLD, _YELLOW, _BOLD, None)
        elif mode == 'resumee':
            colors = (_BOLD, _GREEN, _BOLD, None)
        line = self.build_line(words, colors=colors)
        text = text + line

        # empty line
        line = '\n'
        text = text + line

        # column headers
        words = ['', 'Total', 'Done', 'Remaining']
        colors = (_BOLD,) * len(words)
        line = self.build_line(words, colors=colors)
        if self.print_colors:
            line = _BOLD + line + _ENDC
        text = text + line

        _colors = (_BOLD, None, None, None)

        # time
        words[0] = 'Time:'
        words[1] = utils.time_string(utils.time_total(self))
        words[2] = utils.time_string(utils.time_done(self))
        words[3] = utils.time_string(utils.time_todo(self))
        if mode == 'run':
            colors = (_BOLD, None, None, _YELLOW)
        elif mode == 'resumee':
            colors = (_BOLD, _GREEN, None, None)
        line = self.build_line(words, colors=colors)
        text = text + line

        # count line
        words[0] = 'Loops:'
        words[1] = utils.count_string(self.total_count)
        words[2] = utils.count_string(np.floor(self.count))
        words[3] = utils.count_string(np.ceil(self.total_count - self.count))
        line = self.build_line(words, colors=_colors)
        text = text + line

        # fraction line
        words[0] = 'Fraction:'
        words[1] = '100 %'
        words[2] = string_utils.percentage_string(utils.fraction_done(self))
        words[3] = string_utils.percentage_string(utils.fraction_todo(self))
        line = self.build_line(words, colors=_colors)
        text = text + line

        # bar
        bar_width = sum(_col_width[1:])
        fraction_done = utils.fraction_done(self)
        delim_color = _BOLD
        if mode == 'run':
            fillcolor = ''
        elif mode == 'resumee':
            fillcolor = _GREEN

        bar = string_utils.progress_bar(
                fraction_done, bar_width, fillcolor=fillcolor,
                delim_color=delim_color, use_color=self.print_colors,
                )
        text = text + ' ' * _col_width[0] + bar + '\n'

        # usermessage
        if usermessage is not None:
            python_version = sys.version_info[0]    # (int)
            if python_version < 3:
                usermessage = unicode(usermessage)
            line = '\n' + usermessage + '\n'
            text = text + line

        return text

    def update_screen(self, text):
        self.delete_lines(self.Nlines_last_message)
        self.Nlines_last_message = text.count('\n') + 1
        self._print(text)

    def delete_lines(self, N=0):
        N = int(N)
        if N == 0:
            return self
        if N < 0:
            raise ValueError()

        if self.file is None:
            deleter = N * (_UPWARD + _CLEARLINE) + _UPWARD
            self._print(deleter)
        elif not os.path.isfile(self.file):
            pass
        else:
            with open(self.file, 'r') as fid:
                lines = fid.readlines()
            with open(self.file, 'w') as fid:
                for line in lines[:-N]:
                    fid.write(line)

        return self

    def _initialize_file(self):
        if self.file is None:
            return self

        if self.file_initialized:
            return self

        if os.path.isfile(self.file):
            self.file_initialized = True
            return self

        dirname = os.path.dirname(self.file)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        if not os.path.isfile(self.file):
            with open(self.file, 'w'):
                pass
        self.file_initialized = True

        return self

    def _print(self, text):
        if self.file is None:
            _builtin_print(text)
        else:
            self._initialize_file()
            with open(self.file, 'a') as fid:
                # fid.write((text + '\n').encode('utf-8'))
                fid.write(text + '\n')
        return self

class PerformanceInfo(Chronometer):
    """Alias to Chronometer."""
    pass
