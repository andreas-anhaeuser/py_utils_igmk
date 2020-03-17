import datetime as dt

class MonthTimeDelta(dt.timedelta):
    def __init__(self, number):
        assert isinstance(number, int)
        self.number = number

class YearTimeDelta(dt.timedelta):
    def __init__(self, number):
        assert isinstance(number, int)
        self.number = number

class EnhancedTimeDelta(dt.timedelta):
    def __init__(
            self, days=0, seconds=0, microseconds=0, minutes=0, hours=0,
            weeks=0, years=0, months=0,
            ):
        raise NotImplementedError()
        if years == 0 and months == 0:
            return super().__init__(
                    days=days, seconds=seconds, microseconds=microseconds,
                    minutes=minutes, hours=hours, weeks=weeks,
                    )


####################################################
# test & debug                                     #
####################################################
if __name__ == '__main__':
    print('gg')
