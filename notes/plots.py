import numpy as np
import matplotlib.pyplot as plt

sci = [3, 4, 5, 6, 7, 7.5, 8.0, 9.0, 10.0, 12.0]
y = [49492,46520,43994,42815,41975,41358,41642,40843,40735,41061]
z = [73.8,62.8,64.9,66.7,68.6,71.5,75.9,80.0,82.8,86.2]

import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import datetime
#
# x = [
#     datetime.datetime(2011, 1, 4, 0, 0),
#     datetime.datetime(2011, 1, 5, 0, 0),
#     datetime.datetime(2011, 1, 6, 0, 0)
# ]
x = date2num(x)

#y =  [3.0, 4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0]
z = [49492,46520,43994,42815,41975,41358,41642,40843,40735,41061]
k = [73.8,62.8,64.9,66.7,68.6,71.5,75.9,80.0,82.8,86.2]

ax = plt.subplot(111)
ax.bar(x-0.2, z, width=0.2, color='b', align='center')
ax.bar(x, z, width=0.2, color='g', align='center')
ax.bar(x+0.2, k, width=0.2, color='r', align='center')
ax.xaxis_date()

plt.show()