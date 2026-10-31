"""Reserved source interface. No optical heating is implemented in MR3."""


def heat_source(t, R, T):
    """Power entering the noncondensable gas, W; identically zero here."""
    return 0.0
