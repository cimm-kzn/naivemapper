@julius
def spam(year):
    return None


def julius(spam):
    def wrapper(year):
        if year<0:
            year=-year-1
            spam(year)
        else:
            spam(year)
    return wrapper
