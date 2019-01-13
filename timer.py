#!/usr/bin/python
"""A type that measures elapsed time."""

# standard modules
import datetime as dt

# local modules
if __name__ == '__main__':
    from string_utils import human_format
else:
    from .string_utils import human_format

class Timer(object):
    """Initialize stopped timer.

        Parameters
        ----------
        <none>

        Methods
        -------
        start() : (re-)start counting
        stop() : stop counting
        get() : show elapsed time
        reset() : set to 0. and stop

        Methods may be concatenated:
        >>> timer = Timer().start()
        >>> timer.stop().show().reset().start()

        Examples
        --------
        >>> from time import sleep

        >>> # initialize and start timer
        >>> timer = Timer().start()
        >>> sleep(0.3)

        >>> # show status while running
        >>> print(timer)
        >>> sleep(0.7)

        >>> # show status after stopping
        >>> timer.stop()
        >>> print(timer)

        >>> # show status of stopped timer later
        >>> sleep(0.5)
        >>> print(timer)
    """
    def __init__(self):
        """Return self.
            
            A Timer has two states:
            - running
            - stopped
            When stopped, its time count (field `elapsed`) may be non-zero if
            it has been running previously.

            Fields
            ------
            elapsed : float
                (s) elapsed time while running
            started : None of dt.datetime
                time when it was last started; None if currently stopped.
        """
        self.reset()

    def __str__(self):
        """Return elapsed time as human readable string."""
        secs = self.get()
        return str(human_format(secs)) + 's'

    def __repr__(self):
        """Return a str."""
        elapsed = self.get()
        if self.is_running():
            state = 'running'
        else:
            state = 'stopped'
        return '%s Timer at %f seconds' % (state, elapsed)

    def reset(self):
        """Return a stopped timer with zero elapsed time."""
        self.elapsed = 0.
        self.started = None
        return self

    def start(self, ignore_running=False):
        """Start the timer and return self.
            
            Parameters
            ----------
            ignore_running : bool, optional
                (default: False) If False, an error is thrown, if the timer is
                already running. If True, an already running Timer is left
                unchanged.

            Throws
            ------
            AlreadyRunningError
        """
        if self.is_running():
            if ignore_running:
                return self
            else:
                message = 'Attempt to start already running Timer.'
                raise AlreadyRunningError(message)

        # regular case
        self.started = dt.datetime.now()
        return self

    def stop(self, ignore_stopped=False):
        """Stop the timer and return self.
            
            Parameters
            ----------
            ignore_stopped : bool, optional
                (default: False) If False, an error is thrown, if the timer is
                already stopped. If True, an already stopped Timer is left
                unchanged.

            Throws
            ------
            AlreadyStoppedError
        """
        if self.is_stopped():
            if ignore_stopped:
                return self
            else:
                message = 'Attempt to stop already stopped Timer.'
                raise AlreadyStoppedError(message)

        # regular case
        now = dt.datetime.now()
        diff = (now - self.started).total_seconds()
        self.elapsed += diff
        self.started = None
        return self

    def is_running(self):
        """Return a bool."""
        return not self.is_stopped()

    def is_stopped(self):
        """Return a bool."""
        return self.started is None

    def get(self):
        """Return elapsed time in seconds as float."""
        if self.is_running():
            self.stop()
            self.start()

        return self.elapsed

    def show(self):
        """Print elapsed time and return self."""
        print(self.get())
        return self

    def set(self, elapsed):
        """Set elapsed time and return self.

            Running/stopping state remains unchanged.

            Parameters
            ----------
            elapsed : float or datetime.timedelta
        """
        if isinstance(elapsed, dt.timedelta):
            elapsed = elapsed.total_seconds()

        self.elapsed = float(elapsed)
        return self


class TimerError(Exception):
    pass

class AlreadyRunningError(TimerError):
    pass

class AlreadyStoppedError(TimerError):
    pass

################################################################
# testing                                                      #
################################################################
if __name__ == '__main__':
    timer = Timer()
    timer.start().show()
    timer.set(2).show()
    timer.stop().show().reset().show()