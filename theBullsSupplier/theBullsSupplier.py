import time
import sys
import random
import math
import numpy
import matplotlib.pyplot as pyplot

sys.setrecursionlimit(10000)

#MANDATORY PARAMETERS

#             Begin Date     End Date      Save Plot As
dateList = [['2007-01-01', '2007-12-31', 'strat2007.png'],
            ['2008-01-01', '2008-12-31', 'strat2008.png'],
            ['2009-01-01', '2009-12-31', 'strat2009.png'],
            ['2010-01-01', '2010-12-31', 'strat2010.png'],
            ['2011-01-01', '2011-12-31', 'strat2011.png'],
            ['2012-01-01', '2012-12-31', 'strat2013.png'],
            ['2013-01-01', '2013-12-31', 'strat2013.png']]

maxHoldPeriod = 10 #Do not have to use must be filled out
randomizeMHP = False
numTrials = 10000
trailingPeriods = 200



#PORTFOLIO OPTIONS

initialAmount = 5000
flatRate = 3.95
slippage = .005


#PLOT OPTIONS

plot = True

plotPopulation = True
plotSP500 = True
plotMiddle = True
middlePercent = .25

savePlot = True
showPlot = False


#TICKER OPTIONS

fromFile = True
fileName = 'SP500.txt'
userDefined = False
userList = [' ']


#USER GENERATED PARAMETERS

highlowPeriods = 10
longMA = 200
shortMA = 50


#BUY AND SELL SIGNALS
#Note: Input is a list of dictinaries that contain Open, High, Low, and Close
#data sorted newest first. This function must return a buy/sell price if
#signals are triggered otherwise they must return False.

def buySignal(histData):

    lowoverLongMA = histData[0]['Low'] > sMA(histData[1:], longMA)
    lowoverShortMA = histData[0]['Low'] > sMA(histData[1:], shortMA)
    newLow = histData[0]['Low'] < histLow(histData[1:], highlowPeriods)

    buyPrice = histLow(histData[1:], highlowPeriods)
    if buyPrice > histData[0]['Open']:
        buyPrice = histData[0]['Open']

    if lowoverShortMA and lowoverLongMA and newLow:
        return buyPrice
    else:
        return False


def sellSignal(histData, daysHolding = sys.maxint, maxHoldPeriod = sys.maxint):

    #price crosses ma
    if histData[0]['Low'] < sMA(histData[1:], shortMA):

        if histData[0]['Open'] < sMA(histData[1:], shortMA):
            return histData[0]['Open']
        else:
            return round(sMA(histData[1:], shortMA), 2)

    #new high
    if histHigh(histData[1:], highlowPeriods) < histData[0]['High']:

        if histHigh(histData[1:], highlowPeriods) < histData[0]['Open']:
            return histData[0]['Open']
        else:
            return round(histHigh(histData[1:], highlowPeriods), 2)

    if daysHolding == maxHoldPeriod:
        return histData[0]['Close']
    return False




###############################################################################
##                                                                           ##
##                        BEGINNING OF PROGRAM                               ##
##                                                                           ##
###############################################################################


def portfolio(trade, bank, flatRate):

    #Correct buy price and sell price for slippage
    purchasePrice = trade[2]*(1 + random.uniform(-slippage,slippage))
    sellPrice = trade[4]*(1 + random.uniform(-slippage,slippage))
    purchaseQuant = int(bank/purchasePrice)

    #Subtract comission
    bank = bank - flatRate
    #Subtract cost basis
    bank = bank - (purchaseQuant * purchasePrice)
    #Add market value at sell
    bank = bank + (purchaseQuant * sellPrice)
    #Subtract comission
    bank = bank - flatRate

    return bank

def importData(ticker, beginDate, endDate, trailingPeriods):
    beginDayFound = False
    endDayFound = False
    dataList = []

    f = open("Tickers/" + ticker + ".txt", "r")

    #Skip header line
    line = f.readline().split()
    line = f.readline().split()

    while line:

        #Start collecting on end date
        if int(line[0][:4]+line[0][5:7]+line[0][8:11]) <= int(endDate[:4]+endDate[5:7]+endDate[8:11]):
            endDayFound = True

            #If begin date found collect trailing period data
            if int(line[0][:4]+line[0][5:7]+line[0][8:11]) < int(beginDate[:4]+beginDate[5:7]+beginDate[8:11]):
                beginDayFound = True

                i = 0
                while i < trailingPeriods:

                    #If there is not enough data to support trailing periods return nothing
                    try:
                        tempDict= {}
                        tempDict.update({'Date':line[0]})
                        tempDict.update({'Open':float(line[1])})
                        tempDict.update({'High':float(line[2])})
                        tempDict.update({'Low':float(line[3])})
                        tempDict.update({'Close':float(line[4])})
                        tempDict.update({'Volume':float(line[5])})
                        tempDict.update({'AdjClose':float(line[6])})
                        dataList.append(tempDict)

                    except IndexError:
                        beginDayFound = False

                    line = f.readline().split()
                    i = i + 1
                break

            #If end date found and begin date not found collect data
            else:
                tempDict= {}
                tempDict.update({'Date':line[0]})
                tempDict.update({'Open':float(line[1])})
                tempDict.update({'High':float(line[2])})
                tempDict.update({'Low':float(line[3])})
                tempDict.update({'Close':float(line[4])})
                tempDict.update({'Volume':float(line[5])})
                tempDict.update({'AdjClose':float(line[6])})
                dataList.append(tempDict)

        line = f.readline().split()


    if beginDayFound and endDayFound:
        return dataList
    else:
        return []

#Returns a simple moving average
def sMA(historicalData, periods):
    numerator = 0

    for entry in historicalData[:periods]:
        numerator = numerator + entry['Close']

    return numerator/periods

#Returns the historical low
def histLow(historicalData, periods):
    low = sys.maxint

    for entry in historicalData[:periods]:
        if entry['Low'] < low:
            low = entry['Low']

    return low

#Returns the historical high
def histHigh(historicalData, periods):
    high = 0

    for entry in historicalData[:periods]:
        if entry['High'] > high:
            high = entry['High']

    return high

#This function will collect all possible trades for a single ticker
#and return a list containing lists of all those trades.
#Data will look like this:
#[ticker, buyDate, buyPrice, sellDate, sellPrice, 1+percent gain]
def collectTrades(ticker, historicalData, beginDate):

    tradesList = []

    #Find index of the begin date
    i = 0
    while i < len(historicalData):
        if int(historicalData[i]['Date'][:4]+historicalData[i]['Date'][5:7]+historicalData[i]['Date'][8:11]) < int(beginDate[:4]+beginDate[5:7]+beginDate[8:11]):
            beginDateIndex = i - 1
            break
        i= i + 1

    #Cycle through days
    i = beginDateIndex
    while i >= 0:

        #Check for buy signal
        buyPrice = buySignal(historicalData[i:])
        staticMHP = maxHoldPeriod

        #This will randomly assign a maximum holding period between
        #[0, maxHoldPeriod] and is often used when calculating the
        #returns of base cases
        if randomizeMHP:
            staticMHP = random.randint(1,maxHoldPeriod)

        buyDate = historicalData[i]['Date']

        #if buy signal triggered
        if buyPrice:

            #if trigger occurs on the last day return that days close
            if i == 0:
                tradesList.append([ticker, buyDate, buyPrice, historicalData[i]['Date'], historicalData[i]['Close'], round(historicalData[i]['Close']/buyPrice, 3), 0])
                break

            #oherwise move to next day and cycle through days until sell signal
            #is triggered
            else:
                j = i - 1
                while j >= 0:

                    #check for sell signal
                    sellPrice = sellSignal(historicalData[j:], i - j, staticMHP)
                    sellDate = historicalData[j]['Date']

                    #if sell signal triggered
                    if sellPrice:
                        tradesList.append([ticker, buyDate, buyPrice, sellDate, sellPrice, round(sellPrice/buyPrice, 3), i - j])
                        break

                    #if this position is still open on the last day of analysis
                    #return that days close
                    if not(sellPrice) and j == 0:
                        sellPrice = historicalData[j]['Close']
                        tradesList.append([ticker, buyDate, buyPrice, sellDate, sellPrice, round(sellPrice/buyPrice, 3), i - j])
                        break

                    #calculate next sell signal
                    j = j - 1

        #calculate next buy signal
        i = i - 1

    return tradesList

#Create a trade list for multiple tickers
#Will return a list of lists similar to 'collectTrades'
def universalTradeList(beginDate, endDate, fromFile, userDefined, fileName, userList):

    #import tickers from file
    if fromFile:
        with open(fileName, 'r') as f:
            tickerList = [line.strip() for line in f]

    #user user defined tickers
    if userDefined:
        tickerList = userList

    #collect trades for each ticker and append to completeTrades list
    completeTrades = []
    for ticker in tickerList:

        #collect data for individual tickers
        try:
            historicalData = importData(ticker, beginDate, endDate, trailingPeriods)
        except:
            historicalData = []

        #if there is sufficient data collect trades
        if historicalData:
            for entry in collectTrades(ticker, historicalData, beginDate):
                completeTrades.append(entry)

    return completeTrades

#Swap elements in list
def swap(L, i, j):
    temp = L[i]
    L[i] = L[j]
    L[j] = temp


#Partitions a list of dates for quick sort
def partition(L, first, last, sortColumn):
    # We pick the element L[first] as the "pivot" around which we partition the list
    p = first

    # We process the rest of the elements, one-by-one, in left-to-right order
    for current in range(p+1, last+1):

        #Turn date into int
        currentDate = ''
        for char in L[current][sortColumn]:
            if char != '-':
                currentDate = currentDate + char
        currentDate = int(currentDate)

        #Turn date into int
        pDate= ''
        for char in L[p][sortColumn]:
            if char != '-':
                pDate = pDate + char
        pDate = int(pDate)

        # If L[current] is smaller than the pivot, it needs to move into the first block,
        # to the left of the pivot.
        if currentDate < pDate:
            swap(L, current, p+1)
            swap(L, p, p+1)
            p = p + 1

    return p

#partions a list of ints for quick sort
def partitionInts(L, first, last, sortColumn):
    # We pick the element L[first] as the "pivot" around which we partition the list
    p = first

    # We process the rest of the elements, one-by-one, in left-to-right order
    for current in range(p+1, last+1):

        # If L[current] is smaller than the pivot, it needs to move into the first block,
        # to the left of the pivot.
        currentInt = L[current][sortColumn]
        pInt = L[p][sortColumn]

        if currentInt < pInt:
            swap(L, current, p+1)
            swap(L, p, p+1)
            p = p + 1

    return p


#Recursive quick sort for dates
def generalQuickSort(L, first, last, sortColumn):
    # Base case: if first == last, then there is only one element in the
    # slice that needs sorting. So there is nothing to do.

    # Recursive case: if there are 2 or more elements in the slice L[first:last+1]
    if first < last:
        # Divide step: partition returns an index p such that
        # first <= p <= last and everthing in L[first:p] is <= L[p]
        # and everything in L[p+1:last+1] is >= L[p]
        p = partition(L, first, last,sortColumn)

        # Conquer step
        generalQuickSort(L, first, p-1, sortColumn)
        generalQuickSort(L, p+1, last, sortColumn)

        # Combine step: there is nothing left to do!
    return L

#Recursive quick sort for ints
def generalQuickSortInts(L, first, last, sortColumn):

    # Base case: if first == last, then there is only one element in the
    # slice that needs sorting. So there is nothing to do.

    # Recursive case: if there are 2 or more elements in the slice L[first:last+1]
    if first < last:
        # Divide step: partition returns an index p such that
        # first <= p <= last and everthing in L[first:p] is <= L[p]
        # and everything in L[p+1:last+1] is >= L[p]
        p = partitionInts(L, first, last,sortColumn)

        # Conquer step
        generalQuickSortInts(L, first, p-1, sortColumn)
        generalQuickSortInts(L, p+1, last, sortColumn)

        # Combine step: there is nothing left to do!

    return L


#Main quicksort functino for dates
def quickSort(L, sortColumn):
    return generalQuickSort(L, 0, len(L)-1, sortColumn)

#Merge sort a list of lists that contain quick-sorted integers
#This mitigates the maximum recursion depth errors
def mergeQuickSortInts(L, sortColumn):
    mergedList = []
    buckets = len(L)

    while L:
        minIndex = 0
        minValue = sys.maxint

        #Find the bucket with the lowest value
        for i in range(0, buckets):
            if L[i] and L[i][0][sortColumn] < minValue:
                minValue = L[i][0][sortColumn]
                minIndex = i

        #Append to mergedList and delete from bucket
        if minValue != sys.maxint:
            mergedList.append(L[minIndex][0])
            del L[minIndex][0]
        else:
            break

    return mergedList

#Main quicksort function for ints
def quickSortInts(L, sortColumn):
    buckets = len(L)/1000
    seperatedList = []

    #Seperate into sorted buckets
    i = 0
    while i < buckets:
        temp = L[:1000]
        seperatedList.append(generalQuickSortInts(temp, 0, len(temp) - 1, sortColumn))
        del L[:1000]
        i = i + 1

    #Add any residual buckets to seperatedList
    if len(L) > 1:
        seperatedList.append(generalQuickSortInts(L, 0, len(L) - 1, sortColumn))

    return mergeQuickSortInts(seperatedList,sortColumn)

#This functin will return a dictionary of all of our trades with keywords
#buyDate and values list of possible trades.  Input is a sorted list of trades
#taken from 'universalTradeList'
def createTradeDictionary(L):
    tradeDict = {}
    tradesOnDate = []

    for entry in L:
        if tradesOnDate:

            #If this trade occurred on the same day as those in 'tradesOnDate'
            #append it to 'tradesOnDate'
            if tradesOnDate[0][1] == entry[1]:
                tradesOnDate.append(entry)

            #Otherwise put the list in the dictionary and append the new trade
            #to 'tradesOnDate'
            else:
                tradeDict.update({tradesOnDate[0][1]:tradesOnDate})
                tradesOnDate = []
                tradesOnDate.append(entry)

        #if this is the first trade in L
        else:
            tradesOnDate.append(entry)

    #update dictionary with any residual trades
    if tradesOnDate:
        tradeDict.update({tradesOnDate[0][1]:tradesOnDate})

    return tradeDict

#This function takes our dictionary of trades and uses it to create
#random trade sequences. The output it a list of trades
def createTradeSequence(tradeDict):

    trades = []

    #sort dictionary keys to collect possible trade dates
    tradeDates = []
    for entry in tradeDict.keys():
        tradeDates.append([entry])
    temp = quickSort(tradeDates, 0)
    tradeDates = []
    for entry in temp:
        tradeDates.extend(entry)

    sellDate = 0
    #cycle through possible trade dates
    for entry in tradeDates:

        #make dates ints
        temp = ''
        for char in entry:
            if char != '-':
                temp = temp + char
        modifiedEntry = int(temp)

        #if the current date is after the previous sell date
        if modifiedEntry > sellDate:

            #select a random trade and apend it to 'trades'
            tradeIndex = random.randint(0, len(tradeDict[entry]) - 1)
            trades.append(tradeDict[entry][tradeIndex])

            #format sell date to int and update
            temp = ''
            for char in tradeDict[entry][tradeIndex][3]:
                if char != '-':
                    temp = temp + char
            sellDate = int(temp)

    return trades


#This is the guts of the program and calculates everything we need to know
#about how our strategy performed
def runBacktest(tradeDict, numTrials, plot, testTicker, startDate, endDate):
    returnsList = []
    x = []
    y = []

    #run test for numTrials many times
    aveReturns = 0
    for i in xrange(numTrials):

        #create a trade sequence
        trades = createTradeSequence(tradeDict)

        #collect our portfolio value after each trade
        bank = initialAmount
        for entry in trades:
            bank = portfolio(entry, bank, flatRate)

            #collect our x and y value for plotting
            if plot:
                x.append(entry[3])
                y.append(bank/initialAmount)

        #calculate our trade sequence's return
        returns = bank/initialAmount
        returnsList.append(returns)
        aveReturns += returns

    #calculate our test's average return,sample standard deviation,
    #number of winners, maximum portfolio value, minimum portfolio value

    #average return
    aveReturn = aveReturns/numTrials

    winners = 0
    stdev = 0
    maximum = 0
    minimum = sys.maxint

    for entry in returnsList:
        #number of winners
        if entry >= 1:
            winners += 1
        #maximum return
        if entry > maximum:
            maximum = entry
        #minimumreturn
        if entry < minimum:
            minimum = entry
        stdev +=  (entry - aveReturn)**2
    #standard deviation
    if len(returnsList) > 1:
        stdev = math.sqrt(stdev/(len(returnsList) - 1))

    #calculate number of losers
    losers = numTrials - winners
    #calculate the percent of sequences that won
    percentWinners = float(winners)/numTrials


    #if we're plotting we'll have to reformat x
    if plot:
        #get the days the market was open
        marketDates = []
        test_data = importData(testTicker, beginDate, endDate, 0)
        for day in test_data:
            marketDates.append(day['Date'])

        #put dates in a dictionary with key Date and value x'th trading day of the year
        dateDictionary = {}

        for (i, date) in enumerate(marketDates):
            dateDictionary.update({date: i + 1})

        #reformat x
        tempX = []
        for point in x:
            if point in dateDictionary:
                tempX.append(dateDictionary[point])
        x = tempX

    return [aveReturn, stdev, maximum, minimum, winners, losers, percentWinners, x, y]


#This function will plot the returns of the SP500, you must have
#SPY data downloaded.  Function returns the minimum and maximum y values
#needed for plotting.
def plotSP500(beginDate, endDate):

    x = []
    y = []
    maxY = 0
    minY = sys.maxint

    #import daily SPY data
    dailyData = importData('SPY', beginDate, endDate, 0)
    firstPrice = dailyData[len(dailyData) - 1]['Close']
    profit = 1

    #collect y
    i = len(dailyData) - 1
    while i >= 0:
        point = (dailyData[i]['Close']/firstPrice - 1) * 100
        y.append((round(point, 3)))

        #find min and max y values
        if point > maxY:
            maxY = point
        if point < minY:
            minY = point
        i = i - 1


    #create x
    i = 1
    while i <= len(y):
        x.append(i)
        i = i + 1

    #plot the data
    pyplot.plot(x,y, color = 'red', linewidth = 1.5)

    return(minY, maxY)

#This function takes in a list of [x,y] values and returns
#the middle %
def returnMiddle(pointList, middle):

    #collect all y values
    tempPList = []
    for entry in pointList:
        tempPList.append(entry[1])
    pointList = tempPList

    #find the upper and lower bounds of our middle %
    lowerBound = 0
    upperBound = 0
    if len(pointList) > 1:
        lowerBound = len(pointList)*middle*.5
        upperBound = len(pointList) - lowerBound

    #middle bounds adjusted for index
    lowerBound = int(math.floor(lowerBound))
    upperBound = int(math.ceil(upperBound))

    #collect all data that falls into our range
    listofMiddles = []
    for i in range(lowerBound, upperBound):
        listofMiddles.append(pointList[i])

    return  listofMiddles

#This functino will plot the middle middle% of portfolio values
#on each day a trade is made. The function will return the minimum
#and maximum y values needed for plotting
def plotMiddle(x,y, middle):
    maxY = 0
    minY = sys.maxint

    #create a sorted list that contains lists of points [x,y]
    combinedXY = []
    i = 0
    while i < len(x):
        combinedXY.append([x[i],y[i]])
        i = i + 1
    combinedXY = quickSortInts(combinedXY, 0)

    pointDict = {}
    pointList = []
    currentPoint = ''

    #create a dictionary with Key x and value list of y's that fall on x
    for entry in combinedXY:
        [x,y] = entry

        #if first entry
        if not(pointList):
            pointList.append([x,y])
            currentPoint = x
        else:

            #x matches those in 'pointList' append x to pointList
            if x == currentPoint:
                pointList.append([x,y])

            #otherwise find pointList middle range restart pointList
            else:

                #find midle% of points
                pointList = quickSortInts(pointList, 1)
                pointList = returnMiddle(pointList, middle)

                #insert into pointDict
                pointDict.update({currentPoint: pointList})

                #clear point list and update with new x
                pointList = []
                pointList.append([x,y])
                currentPoint = x

    #repeat for any residual points
    if pointList:

        #insert middle points into dictionary
        pointList = quickSortInts(pointList, 1)
        pointList = returnMiddle(pointList, middle)

        #insert into pointDict
        pointDict.update({currentPoint: pointList})


    #create middle x and y
    middleX = []
    middleY = []
    for xPoint in pointDict:
        for yPoint in pointDict[xPoint]:
            middleX.append(xPoint)
            middleY.append(yPoint)

            #find max and min y
            if yPoint < minY:
                minY = yPoint
            if yPoint > maxY:
                maxY = yPoint

    #draw plot
    pyplot.scatter(middleX,middleY, color = 'blue')

    return [minY, maxY]


#Main plotting function
def plot(x, y, plotFileName, showPlot, savePlot):

    #find x extremes
    maxX = 0
    minX = 0
    for entry in x:
        if entry > maxX:
            maxX = entry

    #plot x axis
    xAxisY = [0]*int(math.ceil(maxX)+1)
    xAxisX = []
    for i in range(0, len(xAxisY)):
        xAxisX.append(i)
    pyplot.plot(xAxisX,xAxisY, color = 'black')

    #reduce y by one
    temp = []
    for point in y:
        temp.append((point - 1) * 100)
    y = temp

    #y extremes
    maxY = 0
    minY = sys.maxint

    #plot population
    if plot and plotPopulation:
        for entry in y:
            if entry > maxY:
                maxY = entry
            if entry < minY:
                minY = entry
        pyplot.scatter(x,y, color = '.80')

    #plot SP500
    if plot and plotSP500:
        [smallY, bigY] = plotSP500(beginDate, endDate)
        if smallY < minY:
            minY = smallY - 5
        if bigY > maxY:
            maxY = bigY + 5

    #plot middle
    if plot and plotMiddle:
        [smallY, bigY] = plotMiddle(x,y, 1 - middlePercent)
        if smallY < minY:
            minY = smallY - 5
        if bigY > maxY:
            maxY = bigY + 5

    #adjust axis scale
    pyplot.axis([minX, maxX, minY, maxY])

    #plot lables
    pyplot.xlabel('Day')
    pyplot.ylabel('Portfolio Value (%)')
    pyplot.title(beginDate + ' - '+ endDate)

    #save plot
    if plot and savePlot:
        pyplot.savefig(plotFileName)

    #view plot
    if plot and showPlot:
        pyplot.show()
    elif plot:
        pyplot.close()

#This is where we put everything together
for entry in dateList:
    start = time.time()
    beginDate = entry[0]
    endDate = entry[1]

    #create tradeList for all tickers in our universe
    List = universalTradeList(beginDate, endDate, fromFile, userDefined, fileName, userList)
    #sort the list by buyDate
    List = quickSort(List, 1)

    #for collecting list of days market was open
    testTicker = None
    if plot:
        testTicker = 'GE'

    #create a trade dictionary
    tradeDict = createTradeDictionary(List)
    print entry, '\n'

    #run our backtest
    output = runBacktest(tradeDict, numTrials, plot, testTicker, beginDate, endDate)
    #[aveReturn, stdev, maximum, minimum, winners, losers, percentWinners, x, y] = runBacktest(tradeDict, numTrials, plot, testTicker, startDate, endDate)

    #print our data
    print 'Average return', output[0]
    print 'Sample STDEV', output[1]
    print 'Maximum return', output[2]
    print 'Minimum return', output[3]
    print 'Winners', output[4]
    print 'Losers', output[5]
    print 'Percent winners', int(output[6]*100), '\n'
    end = time.time()
    print 'Completed in', end - start, '\n'

    #plot our returns
    if True:

        x = output[7]
        y = output[8]
        plot(x,y, entry[2], showPlot, savePlot)

