import inspect

import scipy.signal._savitzky_golay as m

src = inspect.getsource(m.savgol_filter)
i = src.find("fit to the")
print(src[i:i + 900])
