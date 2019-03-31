import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import datetime as dt

from operator import itemgetter

# pd.set_option('display.height', 500)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 300)
np.set_printoptions(threshold=np.nan)

data = pd.read_table("M:\\Lab\\Tyson\\Syne1b master folder\\mutant day-night data\\P22219 s1b-g d196 day-night assay"
                     "\\P22219 s1bg d196 day-night.XLS")

## get subset of data starting at 11pm of night 4 and ending at 10am day 7
# data = data.loc[data["end"] <= 237060]
# data = data.loc[data["end"] >= 24660]
data["sttime"] = pd.to_datetime(data["sttime"])

## convert sttime column to actual time unites in python
#data["sttime"] = data["sttime"].apply(lambda x: dt.datetime.strptime(x, '%H:%M:%S'))
data[["sttime"]].info()
#data["sttime"] = pd.Timestamp(data["sttime"])
print(data["sttime"])

data["sttime"] = data["sttime"].apply(lambda x: pd.Timestamp.round(x, 'min'))
data[["sttime"]].info()
print(data["sttime"])
data["sttime"] = data["sttime"].apply(lambda x: x.time())
data[["sttime"]].info()
print(data["sttime"])
##define a function to round time integers to the nearest minute and call this on sttime to clean up zebrabox times
def roundTime(timeInput=None, dateDelta=dt.timedelta(minutes=1)):
    """Round a datetime object to a multiple of a timedelta
    timeInput : datetime.datetime object, default now.
    dateDelta : timedelta object, we round to a multiple of this, default 1 minute.
    """
    roundTo = dateDelta.total_seconds()

    if timeInput == None : timeInput = dt.datetime.now()
    seconds = (timeInput - timeInput.min).seconds
    # // is a floor division, not a comment on following line:
    rounding = (seconds+roundTo/2) // roundTo*roundTo
    return timeInput + dt.timedelta(0, rounding-seconds)


#rounded_time = data['sttime'].apply(roundTime)
#print(rounded_time)

## This method allows python to generate a list based on lines of a txt file to not manually type them all in
with open("M:\\Lab\\Tyson\\Syne1b master folder\\mutant day-night data\\"
          "P22219 s1b-g d196 day-night assay\\wt_wells.txt", "r") as WW:
    lines = WW.readlines()
    WT_wells = [l.strip() for l in lines if l.strip()]

with open("M:\\Lab\\Tyson\\Syne1b master folder\\mutant day-night data\\"
          "P22219 s1b-g d196 day-night assay\\het_wells.txt", "r") as HW:
    lines = HW.readlines()
    Het_wells = [l.strip() for l in lines if l.strip()]

# with open("M:\\Lab\\Tyson\\Syne1b master folder\\mutant day-night data\\"
#                           "P22219 s1b-g d196 day-night assay\\mutant_wells.txt", "r") as MW:
#   lines = MW.readlines()
# Mut_wells = [l.strip() for l in lines if l.strip()]


# print (WT_wells)
# print(Mut_wells)
# print(Het_wells)
# data_WT = data[data["animal"].isin(WT_wells)]
# data_Mut= data[data["animal"].isin(Mut_wells)]
# data_Het = data[data["animal"].isin(Het_wells)]

# timesplit = data_WT['sttime'].apply(lambda x: x.split(":"))


# print(data_WT)
# print(data_Mut)
# print(data_Het)
##get data divided out into each day and night these numbers are calculated based on the end time of the data in sec
WT_night4 = data_WT.loc[data_WT["end"] <= 63900]
WT_day5 = data_WT.loc[data_WT["end"] >= 63960]
WT_day5 = WT_day5.loc[WT_day5["end"] <= 114300]
WT_night5 = data_WT.loc[data_WT["end"] >= 114360]
WT_night5 = WT_night5.loc[WT_night5["end"] <= 150300]
WT_day6 = data_WT.loc[data_WT["end"] >= 150360]
WT_day6 = WT_day6.loc[WT_day6["end"] <= 200700]
WT_night6 = data_WT.loc[data_WT["end"] >= 200760]
WT_night6 = WT_night6.loc[WT_night6["end"] <= 236700]

lg_night4 = data_lg.loc[data_lg["end"] <= 63900]
lg_day5 = data_lg.loc[data_lg["end"] >= 63960]
lg_day5 = lg_day5.loc[lg_day5["end"] <= 114300]
lg_night5 = data_lg.loc[data_lg["end"] >= 114360]
lg_night5 = lg_night5.loc[lg_night5["end"] <= 150300]
lg_day6 = data_lg.loc[data_lg["end"] >= 150360]
lg_day6 = lg_day6.loc[lg_day6["end"] <= 200700]
lg_night6 = data_lg.loc[data_lg["end"] >= 200760]
lg_night6 = lg_night6.loc[lg_night6["end"] <= 236700]

sh_night4 = data_sh.loc[data_sh["end"] <= 63900]
sh_day5 = data_sh.loc[data_sh["end"] >= 63960]
sh_day5 = sh_day5.loc[sh_day5["end"] <= 114300]
sh_night5 = data_sh.loc[data_sh["end"] >= 114360]
sh_night5 = sh_night5.loc[sh_night5["end"] <= 150300]
sh_day6 = data_sh.loc[data_sh["end"] >= 150360]
sh_day6 = sh_day6.loc[sh_day6["end"] <= 200700]
sh_night6 = data_sh.loc[data_sh["end"] >= 200760]
sh_night6 = sh_night6.loc[sh_night6["end"] <= 236700]

##get area under curve for each fish
WT_night4_area = WT_night4.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
WT_day5_area = WT_day5.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
WT_night5_area = WT_night5.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
WT_day6_area = WT_day6.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
WT_night6_area = WT_night6.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
lg_night4_area = lg_night4.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
lg_day5_area = lg_day5.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
lg_night5_area = lg_night5.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
lg_day6_area = lg_day6.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
lg_night6_area = lg_night6.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
sh_night4_area = sh_night4.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
sh_day5_area = sh_day5.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
sh_night5_area = sh_night5.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
sh_day6_area = sh_day6.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)
sh_night6_area = sh_night6.groupby('animal', as_index=True)['actinteg'].agg(np.trapz)

print(WT_night4_area)
all_night4 = []
all_night4.extend(WT_night4_area_821)
all_night4.extend(WT_night4_area)
print(all_night4)

##do t test and mann whitney test comparing days and nights of area under curve each group
print("stats for WTvlg night 4")
print(stats.mannwhitneyu(WT_night4_area, lg_night4_area))
print("stats for WTvlg day 5")
print(stats.mannwhitneyu(WT_day5_area, lg_day5_area))
print("stats for WTvlg night 5")
print(stats.mannwhitneyu(WT_night5_area, lg_night5_area))
print("stats for WTvlg day6")
print(stats.mannwhitneyu(WT_day6_area, lg_day6_area))
print("stats for WTvlg night 6")
print(stats.mannwhitneyu(WT_night6_area, lg_night6_area))

print("stats for WTvsh night 4")
print(stats.mannwhitneyu(WT_night4_area, sh_night4_area))
print("stats for WTvsh day 5")
print(stats.mannwhitneyu(WT_day5_area, sh_day5_area))
print("stats for WTvsh night 5")
print(stats.mannwhitneyu(WT_night5_area, sh_night5_area))
print("stats for WTvsh day6")
print(stats.mannwhitneyu(WT_day6_area, sh_day6_area))
print("stats for WTvsh night 6")
print(stats.mannwhitneyu(WT_night6_area, sh_night6_area))

# data taking total movement for each fish during each period


## generate box and whisker plots to show the values for area under the curve, total movement, and ave movement per min


plt.figure(1)
night4_area = [WT_night4_area, lg_night4_area, sh_night4_area]
plt.boxplot(night4_area, whis=[5, 95], labels=["WT", "lg", "sh"], positions=[1, 2, 3], widths=0.5)
plt.xlim(0.25, 3.75)
plt.title("Area Under the Curve Night 4")

plt.figure(2)
day5_area = [WT_day5_area, lg_day5_area, sh_day5_area]
plt.boxplot(day5_area, whis=[5, 95], labels=["WT", "lg", "sh"], positions=[1, 2, 3], widths=0.5)
plt.xlim(0.25, 3.75)
plt.title("Area Under the Curve Day 5")

plt.figure(3)
night5_area = [WT_night5_area, lg_night5_area, sh_night5_area]
plt.boxplot(night5_area, whis=[5, 95], labels=["WT", "lg", "sh"], positions=[1, 2, 3], widths=0.5)
plt.xlim(0.25, 3.75)
plt.title("Area Under the Curve Night 5")

plt.figure(4)
day6_area = [WT_day6_area, lg_day6_area, sh_day6_area]
plt.boxplot(day6_area, whis=[5, 95], labels=["WT", "lg", "sh"], positions=[1, 2, 3], widths=0.5)
plt.xlim(0.25, 3.75)
plt.title("Area Under the Curve Day 6")

plt.figure(5)
night6_area = [WT_night6_area, lg_night6_area, sh_night6_area]
plt.boxplot(night6_area, whis=[5, 95], labels=["WT", "lg", "sh"], positions=[1, 2, 3], widths=0.5)
plt.xlim(0.25, 3.75)
plt.title("Area Under the Curve Night 6")

# plt.show()

# print(WT_day5_act_min)
# print(WT_day5_act)

## calculate the average movement by minute and standard errors for each group starting at midnight
data_WT_midnight = data_WT.loc[data_WT["end"] >= 31560]
data_lg_midnight = data_lg.loc[data_lg["end"] >= 31560]
data_sh_midnight = data_sh.loc[data_sh["end"] >= 31560]
WT_act_by_min = data_WT_midnight.groupby('end', as_index=True)['actinteg'].agg([np.mean, stats.sem])
lg_act_by_min = data_lg_midnight.groupby('end', as_index=True)['actinteg'].agg([np.mean, stats.sem])
sh_act_by_min = data_sh_midnight.groupby('end', as_index=True)['actinteg'].agg([np.mean, stats.sem])
# print(WT_act_by_min)

## generate arrays containing x and y values to plot and arrays with standard errors
x = np.array(range(526, len(WT_act_by_min) + 526))
# print (x)

y = np.around(WT_act_by_min['mean'].values, decimals=2)
y2 = np.around(lg_act_by_min['mean'].values, decimals=2)
y3 = np.around(sh_act_by_min['mean'].values, decimals=2)

# print(y)

error = np.around(WT_act_by_min['sem'].values, decimals=2)
error2 = np.around(lg_act_by_min['sem'].values, decimals=2)
error3 = np.around(sh_act_by_min['sem'].values, decimals=2)

plt.figure(6)
plt.plot(x, y, 'k-', linewidth=0.25, color="b")
plt.fill_between(x, y - error, y + error, facecolor='#add8e6')

# plt.plot(x, y2, 'k-', linewidth=0.25, color="r")
# plt.fill_between(x, y2-error2, y2+error2, facecolor='#f8cdcf')

plt.plot(x, y3, 'k-', linewidth=0.25, color="g")
plt.fill_between(x, y3 - error3, y3 + error3, facecolor='#b1e0c4')

plt.ylim(-10, 750)

plt.show()
