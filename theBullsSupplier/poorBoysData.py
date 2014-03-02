#theBrokeQuant
#Last Modified: 4/29/2013

import os
import urllib
import sys
from datetime import date

StartDay = 1
StartMonth = 'January'
StartYear = 2005

today = date.today()
EndDay = today.day
EndMonth = today.strftime('%B')
EndYear = today.year

TimeInterval = "Day"

adjustPrices = True

fromFile = True
fileName = 'SP500.txt'

isUsersTickers = False
usersTickerList = ['SPY']





def importData(StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, Ticker, TimeInterval):

    #Obtaining index value for months
    OutputList = []
    Months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    beginmonth = 0
    endmonth = 0
    for i in range (0,12):
        if Months[i] == StartMonth:
            beginmonth = i
        if Months[i] == EndMonth:
            endmonth = i
    #first list element is the beginmonth index, second is endmonth
    OutputList.append(beginmonth)
    OutputList.append(endmonth)

    #construct the URL
    URL = "http://ichart.finance.yahoo.com/table.csv?s=" + Ticker + "&a=" + str(OutputList[0]) + "&b=" + str(StartDay) + "&c=" + str(StartYear) + "&d=" + str(OutputList[1]) + "&e=" + str(EndDay) + "&f=" + str(EndYear)+ "&g=" + TimeInterval[0] + "&ignore=.txt"
    f = urllib.urlopen(URL)
    #save to a temp file
    file("temp.txt", "w").write(f.read())

    #import temp contents into list L
    L = []
    File = open('temp.txt','rt')
    FileList = File.readlines()
    for filerow in range (0, len(FileList)):
        L.append(str(FileList[filerow]).split(","))
    File.close()

    #delete temp file
    os.remove('temp.txt')

    return L

#write the list to a tab seperated file
def listToTabSeperatedFile(L, filename):
    f = open(filename, "w+")
    for entry in L:
        string = ""

        i = 0
        while i < len(entry):
            if i == 0:
                string = string + str(entry[i])
            else:
                string = string +"  " + str(entry[i])
            i = i + 1
        if string[0] == 'D':
            f.write(string)
        else:
            f.write(string + '\n')
    f.close()

tickerList = []

if fromFile:
    #open SP500 list and import tickers
    f= open(fileName,"r")
    lines = f.readline()
    while lines:
        tickerList.append(lines.rstrip())
        lines = f.readline()
    f.close()
if isUsersTickers:
    tickerList = usersTickerList


for Ticker in tickerList:

    try:
        #import the data for the ticker
        L = importData(StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, Ticker, TimeInterval)
        #if we're adjusting prices
        if adjustPrices:
            newL = []
            for entry in L:
                if entry[0][0] != 'D':
                    adjClose = float((entry[6].split("/n"))[0])
                    #open = nonadjOpen/nonadjustClose * adjClose
                    adjOpen = str(round((float(entry[1])/float(entry[4]))*adjClose,2))
                    #high = nonadjHigh/nonadjustClose * adjClose
                    adjHigh = str(round((float(entry[2])/float(entry[4]))*adjClose,2))
                    #low = nonadjLow/nonadjustClose * adjClose
                    adjLow =  str(round((float(entry[3])/float(entry[4]))*adjClose,2))
                    tempL = [entry[0], adjOpen, adjHigh, adjLow, adjClose, entry[5], adjClose]
                    newL.append(tempL)
                else:
                    newL.append(entry)

        L = newL


        #save data to file Ticker.txt in the folder Tickers asa tab seperated file
        filename = "Tickers/" + Ticker + ".txt"
        listToTabSeperatedFile(L, filename)
        print Ticker
    except:
        print "Unknown error, skipped", Ticker





