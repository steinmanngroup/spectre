class SpectreValueError(ValueError):
    pass


class SpectrePotentialValueError(SpectreValueError):
    pass


class SpectreExcitedStateValueError(SpectreValueError):
    pass


class SpectreRuntimeError(RuntimeError):
    pass


class SpectreFolderNotFoundError(FileNotFoundError):
    pass

class SpectreLopropFolderNotFoundError(FileNotFoundError):
    pass


class SpectreLopropFileNotFoundError(FileNotFoundError):
    pass

class SpectrePEEXFileNotFoundError(FileNotFoundError):
    pass