from datetime import datetime

from matplotlib import pyplot as plt

import testUtil


dateToLoad = datetime(2000, 1, 2)
timeIndex = -1
pt = [6., 34.]

# in about 2% of the points the results are strange (1 peak divided into 2) e.g.
#dateToLoad = datetime(2000, 1, 1, 18, 0, 0)
#timeIndex = 17
#pt = [7.80206835, 30.39962121]



outNc = 'outputs/schout_1.nc'

p1hsPt, p1tmPt, p1dirmPt, p1dsprPt, p2hsPt, p2tmPt, p2dirmPt, p2dsprPt, p3hsPt, p3tmPt, p3dirmPt, p3dsprPt = testUtil.getPeaksAtPoint(outNc, timeIndex, pt)

tts, drs, spec = testUtil.getSpecAtDate('P-1.sp2d', dateToLoad)

plt.figure()
plt.pcolor(tts, drs, spec)
plt.scatter([p1tmPt, p2tmPt, p3tmPt], [p1dirmPt, p2dirmPt, p3dirmPt], s=[10*p1hsPt**2, 10*p2hsPt**2, 10*p3hsPt**2])
plt.colorbar()

plt.figure()
plt.pcolor(tts**-1, drs, spec)
plt.scatter([p1tmPt**-1, p2tmPt**-1, p3tmPt**-1], [p1dirmPt, p2dirmPt, p3dirmPt], s=[10*p1hsPt**2, 10*p2hsPt**2, 10*p3hsPt**3])
plt.colorbar()

plt.show()


