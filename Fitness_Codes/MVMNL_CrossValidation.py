#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: chengyilyu
"""

import pandas as pd
import numpy as np
import math
from scipy.optimize import minimize
import random
import datetime

def GammaValue(utilityArray):
    global itCnt,prodAisleMatrixTrain,orderMatrixTrain,prodUtiMatrix, orderMatrixTest, prodAisleMatrixTest,globalNoPurchase
    Gamma = 1
    for key in prodAisleMatrixTrain.keys():
        tempgamma = 1
        for product in prodAisleMatrixTrain[key]:
            index = prodUtiMatrix[product]
            if index == -1:
                tempgamma += np.exp(M)
            else:
                tempgamma += np.exp(utilityArray[index])
        Gamma *= tempgamma
    return Gamma

def callbackFMVMNL(utiltiyArray_k):
    global itCnt,prodAisleMatrixTrain,orderMatrixTrain,prodUtiMatrix, orderMatrixTest, prodAisleMatrixTest,globalNoPurchase
    fnValue = neg_loglikeMVMNL(utiltiyArray_k)
    s = 'iter='+str(itCnt)+' LogLikelihood='+str(fnValue)+'\n',
    print(s)
    itCnt += 1
    return

def neg_loglikeMVMNL(utilityArray, *args):
    global itCnt,prodAisleMatrixTrain,orderMatrixTrain,prodUtiMatrix, orderMatrixTest, prodAisleMatrixTest,globalNoPurchase
    Gamma = GammaValue(utilityArray)
    utility0 = math.log(globalNoPurchase / (1 - globalNoPurchase) * (Gamma - 1))
    negLL = 0
    for key in orderMatrixTrain:
        for product in orderMatrixTrain[key]:
            index = prodUtiMatrix[product]
            if index == -1:
                negLL += M
            else:
                negLL += utilityArray[index]
            negLL -= math.log(np.exp(utility0) - 1 + Gamma)
    negLL = -negLL
    return negLL

def predictChoiceProbability(utilityArray):
    global itCnt,prodAisleMatrixTrain,orderMatrixTrain,prodUtiMatrix, orderMatrixTest, prodAisleMatrixTest,globalNoPurchase
    Gamma = 1
    predictedProba = dict()
    for key in prodAisleMatrixTest.keys():
        tempgamma = 1
        for product in prodAisleMatrixTest[key]:
            index = prodUtiMatrix[product]
            if index == -1:
                tempgamma += np.exp(M)
            else:
                tempgamma += np.exp(utilityArray[index])
        Gamma *= tempgamma
    utility0 = math.log(globalNoPurchase / (1 - globalNoPurchase) * (Gamma - 1))
    for orderIDTest in orderMatrixTest.keys():
        tempProb = 1
        for orderProdTest in orderMatrixTest[orderIDTest]:
            index = prodUtiMatrix[orderProdTest]
            if index == -1:
                tempProb *= np.exp(M)
            else:
                tempProb *= np.exp(utilityArray[index])
        predictedProba[orderIDTest] = tempProb/(np.exp(utility0) - 1 + Gamma)
    return predictedProba

#%%
global itCnt,prodAisleMatrixTrain,orderMatrixTrain,prodUtiMatrix, orderMatrixTest, prodAisleMatrixTest,globalNoPurchase
path = "E:\\phd\\3_Multi-item\\Instacart 2017\\"
# path = "D:\\Users\\chly4234\\Desktop\\3_Multi-item\\Instacart 2017\\"
# Product_count(path)

# dfOrder = pd.read_csv(path + "orders_train.csv",header=0)
dfOrderOrigin = pd.read_csv(path + "orders_test.csv",header=0)
dfAislesOrigin = pd.read_csv(path + "aisles.csv",header=0)
dfProductsOrigin = pd.read_csv(path + "ProductsWithCount.csv",header=0)
dfProducts = dfProductsOrigin.copy()
dfProducts = dfProducts[dfProducts["product_count"]>=20] # only consider the products with at least 35 orders
# utility0 = math.log(dfProducts.shape[0]/2)
globalNoPurchase = 0.3
M = -300

tempOrderProdArray = np.array(dfOrderOrigin["product_id"])
tempProArray = np.array(dfProducts["product_id"])
dfOrder = dfOrderOrigin.copy()
for index in range(dfOrderOrigin.shape[0]):
    tempOrderProd = tempOrderProdArray[index]
    if tempOrderProd not in tempProArray:
        dfOrder.drop([index], inplace = True)
        continue

random.seed(datetime.datetime.now())
samplePathNum = 1
foldnum = 5
randomOrder = [num for num in range(len(dfOrder))]
random.shuffle(randomOrder)
foldsize = int(len(dfOrder)/foldnum)

averageMSE, averageChiSquared, averageKLvalue, totalMSE, totalChiSquared, totalKLvalue = 0.0, 0.0, 0.0, 0.0,0.0,0.0
averageRevenueMSE, totalRevenueMSE = 0.0, 0.0
outofsamplenegLL, averageoutofsamplenegLL, totaloutofsamplenegLL = 0.0, 0.0, 0.0
for path in range(samplePathNum):
    for k in range(foldnum):
        obsIndexTrain, obsIndexTest = [], []
        # pos = 0
        for index in range(len(dfOrder)):
            if (np.float(k) * foldsize <= randomOrder[index]) and (randomOrder[index] < np.float((k + 1)) * foldsize):
                obsIndexTest.append(index)
            else:
                obsIndexTrain.append(index)
            # pos += 1
        orderTrain = dfOrder.iloc[obsIndexTrain,:]
        orderTest = dfOrder.iloc[obsIndexTest,:]

        orderMatrixTrain = dict()
        prodAisleMatrixTrain = dict()
        currentAisleIDTrain = []
        currentOrdIDTrain = []
        currentProdTrain= []
        tempOrdIDArrayTrain = np.array(orderTrain["order_id"])
        tempOrderProdArrayTrain = np.array(orderTrain["product_id"])
        tempOrderAisArrayTrain = np.array(orderTrain["aisle_id"])

        for index in range(orderTrain.shape[0]):
            tempOrdID = tempOrdIDArrayTrain[index]
            tempOrderProd = tempOrderProdArrayTrain[index]
            tempOrderAis = tempOrderAisArrayTrain[index]
            if tempOrderProd not in currentProdTrain:
                currentProdTrain.append(tempOrderProd)
            if tempOrdID not in currentOrdIDTrain:
                currentOrdIDTrain.append(tempOrdID)
                orderMatrixTrain[tempOrdID] = []
                orderMatrixTrain[tempOrdID].append(tempOrderProdArrayTrain[index])
            else:
                orderMatrixTrain[tempOrdID].append(tempOrderProdArrayTrain[index])
            if tempOrderAis not in currentAisleIDTrain:
                currentAisleIDTrain.append(tempOrderAis)
                prodAisleMatrixTrain[tempOrderAis] = []
                prodAisleMatrixTrain[tempOrderAis].append(tempOrderProdArrayTrain[index])
            else:
                if tempOrderProdArrayTrain[index] not in prodAisleMatrixTrain[tempOrderAis]:
                    prodAisleMatrixTrain[tempOrderAis].append(tempOrderProdArrayTrain[index])

        itCnt = 1
        utilityStart = 0.01 * np.random.rand(len(currentProdTrain))
        prodUtiMatrix = dict()
        for index in range(len(currentProdTrain)):
            key = currentProdTrain[index]
            prodUtiMatrix[key] = index
        res = minimize(neg_loglikeMVMNL, utilityStart, method="L-BFGS-B", jac=None, hess=None, callback=callbackFMVMNL,options={'disp': True, 'maxiter': 100})
        utilityArray = res.x
        M = np.min(utilityArray)

        tempOrderProdArrayTest = np.array(orderTest["product_id"])
        orderTestClean = orderTest.copy()
        orderTestClean.index = range(len(orderTestClean))
        for index in range(len(orderTest)):
            tempOrderProd = tempOrderProdArrayTest[index]
            if tempOrderProd not in currentProdTrain:
                prodUtiMatrix[tempOrderProd] = -1

        orderMatrixTest = dict()
        prodAisleMatrixTest = dict()
        currentAisleIDTest = []
        currentOrdIDTest = []
        tempOrdIDArrayTest = np.array(orderTestClean["order_id"])
        tempOrderProdArrayTest = np.array(orderTestClean["product_id"])
        tempOrderAisArrayTest = np.array(orderTestClean["aisle_id"])

        for index in range(orderTestClean.shape[0]):
            tempOrdID = tempOrdIDArrayTest[index]
            tempOrderProd = tempOrderProdArrayTest[index]
            tempOrderAis = tempOrderAisArrayTest[index]
            if tempOrdID not in currentOrdIDTest:
                currentOrdIDTest.append(tempOrdID)
                orderMatrixTest[tempOrdID] = []
                orderMatrixTest[tempOrdID].append(tempOrderProdArrayTest[index])
            else:
                orderMatrixTest[tempOrdID].append(tempOrderProdArrayTest[index])
            if tempOrderAis not in currentAisleIDTest:
                currentAisleIDTrain.append(tempOrderAis)
                prodAisleMatrixTest[tempOrderAis] = []
                prodAisleMatrixTest[tempOrderAis].append(tempOrderProdArrayTest[index])
            else:
                if tempOrderProdArrayTest[index] not in prodAisleMatrixTest[tempOrderAis]:
                    prodAisleMatrixTest[tempOrderAis].append(tempOrderProdArrayTest[index])

        GammaTest = 1
        for key in prodAisleMatrixTest.keys():
            tempgamma = 1
            for product in prodAisleMatrixTest[key]:
                index = prodUtiMatrix[product]
                if index == -1:
                    outofsamplenegLL += np.exp(M)
                else:
                    tempgamma += np.exp(utilityArray[index])
            GammaTest *= tempgamma
        utility0Test = math.log(globalNoPurchase / (1 - globalNoPurchase) * (GammaTest - 1))
        for key in orderMatrixTest:
            for product in orderMatrixTest[key]:
                index = prodUtiMatrix[product]
                if index == -1:
                    outofsamplenegLL += M
                else:
                    outofsamplenegLL += utilityArray[index]
                outofsamplenegLL -= math.log(np.exp(utility0Test) - 1 + GammaTest)
        outofsamplenegLL = -outofsamplenegLL
        totaloutofsamplenegLL += outofsamplenegLL
        print("totaloutofsamplenegLL in (path,k):", path, k, outofsamplenegLL)
        print("totaloutofsamplenegLL:", totaloutofsamplenegLL)

        # Prediction using Test Data
        prob1 = predictChoiceProbability(utilityArray)  # dictionary {"orderID","ChoiceProbability")}

        # KL
        KLvalue1 = 0
        for key in prob1.keys():
            pred, actual = prob1[key], 1.0
            KLvalue1 += actual * np.log(actual / pred)
        # if len(orderTest) == 0:
        #     skipNum += 1
        #     continue
        KLvalue = KLvalue1 / len(prob1)
        totalKLvalue += KLvalue
        print("KLvalue in (path,k):", path, k, KLvalue)
        print("totalKLvalue:", totalKLvalue)

        # Chi-Squared
        ChiSquared1 = 0.0
        for key in prob1.keys():
            pred, actual = prob1[key], 1.0
            ChiSquared1 += (np.square(actual - pred) / pred)
        ChiSquared = ChiSquared1 / len(prob1)
        totalChiSquared += ChiSquared

        ## Revenue MSE
        # RevenueMSE1 = 0.0
        # for i in range(len(orderTest)):
        #     tempobsRevenue, temppredRevenue = 0.0, 0.0
        #     pred, actual = np.array(prob1[i]), np.array(prob2[i])
        #     predsetPrice = np.array(orderTest.iloc[i,14:])
        #     tempobsRevenue = sum(predsetPrice * actual)
        #     temppredRevenue = sum(predsetPrice * pred)
        #     RevenueMSE1 += np.square(tempobsRevenue - temppredRevenue)
        # RevenueMSE = RevenueMSE1 / len(orderTest)
        # totalRevenueMSE += RevenueMSE

    # averageRevenueMSEnueMSE = totalRevenueMSE / (foldnum*samplePathNum - skipNum)
    # print("averageRevenueMSEnueMSE:", averageRevenueMSEnueMSE)

averageKLvalue = totalKLvalue / (foldnum * samplePathNum)
print("averageKL:", averageKLvalue)

averageChiSquared = totalChiSquared / (foldnum * samplePathNum)
print("averageChiSquared:", averageChiSquared)

averageoutofsamplenegLL = totaloutofsamplenegLL / (foldnum * samplePathNum)
print("averageoutofsamplenegLL:", averageoutofsamplenegLL)




