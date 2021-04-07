import numpy as np
import socket
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import pickle
import _pickle as cPickle
from event import Event
import tkinter as Tk
from tkinter import filedialog as tkFileDialog
import time
from matplotlib.patches import Rectangle
from colour import Color
import os

import threading
import datetime
from string import ascii_lowercase as alphabet

matplotlib.use("TkAgg")

#Colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
Colors = [u'b', u'g', u'r', u'c', u'm', u'y', u'k']

def KillDisplay(mainPort = 54242, questionPort = 54243, address = "localhost"):
    Question = b"kill#"
    QuestionUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    QuestionAddress = (address, questionPort)
    QuestionUDP.sendto(Question,QuestionAddress)
    QuestionUDP.close()

    Main = b"kill#"
    MainUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    MainAddress = (address, mainPort)
    MainUDP.sendto(Main,MainAddress)
    MainUDP.close()

class NullHandler:
    DisplayType = 'None'
    def __init__(self):
        pass
    def Decrypt(self, Socket, t, eventList):
        pass
    def AddStreamVars(self, Name, Geometry):
        pass
    def DelStreamVars(self, Name):
        pass
    def Update(self, ReloadMap):
        pass

class EventHandler:
    Key = 1
    DisplayType = "Map"
    Name = "Events"
    def __init__(self, Master, Display, Canvas, InfoFrame, OptionsFrame, SharedData):
        self.Ax = Display
        self.Canvas = Canvas
        self.Active = ((self.Key == 1) or (self.DisplayType == "Dot"))
        self.SharedData = SharedData
        self.Reload = SharedData['ReloadCommand']

        self.DisplayCMap = 'binary'
        self.VMax = 2.
        self.StreamsShapes = {}
        self.StreamsMaps = {}

        L = Tk.Label(InfoFrame, text = "Events shown:")
        L.pack(anchor = Tk.W)
        self.EventsLabel = Tk.Label(InfoFrame)
        self.EventsLabel.pack(anchor = Tk.W)

        PolaritiesFrame = Tk.Frame(OptionsFrame)
        PolaritiesFrame.grid(row = 0, column = 0, sticky = Tk.W)

        PolasModes = [("OFF", 0),
                     ("ON", 1)]
        self.DisplayedPolas = {Polarity: True for PolaName, Polarity in PolasModes}
        self.DisplayedPolaritiesVars = []

        PolaritiesButtons = []
        for nPola, PolaParams in enumerate(PolasModes):
            Text, Polarity = PolaParams
            self.DisplayedPolaritiesVars += [Tk.IntVar(master = Master)]
            self.DisplayedPolaritiesVars[-1].set(int(self.DisplayedPolas[Polarity]))
            PolaritiesButtons += [Tk.Checkbutton(PolaritiesFrame, text = Text, variable = self.DisplayedPolaritiesVars[-1], command = lambda n=nPola, P=Polarity: self.SwitchDisplayedPolasOnOff(n,P))]
            PolaritiesButtons[-1].grid(row = 0, column = nPola, sticky = Tk.N)
        
    def Decrypt(self, Socket, t, eventList):
        self.StreamsMaps[Socket][eventList[0][0], eventList[0][1], eventList[1]] = t
    def AddStreamVars(self, Name, Geometry):
        self.StreamsShapes[Name] = Geometry
        self.StreamsMaps[Name] = -10*np.ones(self.StreamsShapes[Name])
    def DelStreamVars(self, Name):
        del self.StreamsMaps[Name]
        del self.StreamsShapes[Name]
    def Update(self, Reload = False):
        CurrentStream = self.SharedData['CurrentStream']
        if self.Active:
            Map = np.zeros(self.StreamsShapes[CurrentStream][:2])
            StreamMap = self.StreamsMaps[CurrentStream]
            StreamTime = self.SharedData['t']
            for Polarity, IsActive in self.DisplayedPolas.items():
                if IsActive:
                    Map = (Map + ((StreamTime - StreamMap[:,:,Polarity]) < self.SharedData['Tau'])) #For color mode

            if Reload:
                self.DisplayImShow = self.Ax.imshow(np.transpose(Map), vmin = 0, vmax = self.VMax, origin = "lower", cmap = self.DisplayCMap)
                self.Ax.set_xlim(-0.5, self.StreamsShapes[CurrentStream][0]-0.5)
                self.Ax.set_ylim(-0.5, self.StreamsShapes[CurrentStream][1]-0.5)
            else:
                self.DisplayImShow.set_data(np.transpose(Map))
            self.EventsLabel['text'] = "{0}".format(int(Map.sum()))
        else:
            if Reload:
                self.EventsLabel['text'] = "0"

    def SwitchDisplayedPolasOnOff(self, nPola, Polarity):
        if self.DisplayedPolaritiesVars[nPola].get():
            self.DisplayedPolas[Polarity] = True
        else:
            self.DisplayedPolas[Polarity] = False
        self.Reload()

class DisparityHandler:
    Key = 3
    DisplayType = "Map"
    Name = "Disparities"
    def __init__(self, Master, Display, Canvas, InfoFrame, OptionsFrame, SharedData):
        self.Ax = Display
        self.Canvas = Canvas
        self.Active = ((self.Key == 1) or (self.DisplayType == "Dot"))
        self.SharedData = SharedData
        self.Reload = SharedData['ReloadCommand']
        
        self.NCOLORS = 10

        self.MinDisparity = 0
        self.MaxDisparity = 50
        self.DisparitiesMaps = {}
        self.StreamsShapes = {}
        self.Compensate = {}

        L = Tk.Label(InfoFrame, text = "Disparities shown:")
        L.grid(row = 0, column = 0, sticky = Tk.NW)
        self.EventsLabel = Tk.Label(InfoFrame)
        self.EventsLabel.grid(row = 0, column = 1, sticky = Tk.NW)
        self.EventsLabel['text'] = "0"

        CBFig = Figure(figsize=(4,0.5), dpi=60)
        CBFig.subplots_adjust(bottom=0.5)
        self.CBCanvas = FigureCanvasTkAgg(CBFig, InfoFrame)
        CBFig.tight_layout()
        self.CBCanvas.get_tk_widget().grid(row = 1, column = 0, columnspan = 2, sticky = Tk.NW)

        self.CBAx = CBFig.add_subplot(111)
        self.DrawColorbar()

        self.MaxDisparityLabel = Tk.Label(OptionsFrame, text = "Max Disparity : {0}".format(self.MaxDisparity))
        self.MaxDisparityLabel.grid(column = 0, row = 0)
        MaxDMinusButton = Tk.Button(OptionsFrame, text = ' - ', command = lambda: self.ChangeMaxD(-1))
        MaxDMinusButton.grid(column = 1, row = 0)
        MaxDPlusButton = Tk.Button(OptionsFrame, text = ' + ', command = lambda: self.ChangeMaxD(+1))
        MaxDPlusButton.grid(column = 2, row = 0)

        self.MinDisparityLabel = Tk.Label(OptionsFrame, text = "Min Disparity : {0}".format(self.MinDisparity))
        self.MinDisparityLabel.grid(column = 0, row = 1)
        MinDMinusButton = Tk.Button(OptionsFrame, text = ' - ', command = lambda: self.ChangeMinD(-1))
        MinDMinusButton.grid(column = 1, row = 1)
        MinDPlusButton = Tk.Button(OptionsFrame, text = ' + ', command = lambda: self.ChangeMinD(+1))
        MinDPlusButton.grid(column = 2, row = 1)

        AutoDisparityButton = Tk.Button(OptionsFrame, text = "Find disparity range", command = self.AutoSet)
        AutoDisparityButton.grid(column = 0, row = 2)

        self.CompensateVar = Tk.IntVar(master = Master)
        CompensateButton = Tk.Checkbutton(OptionsFrame, text = 'Compensate(c)', variable = self.CompensateVar, command = self.SwitchCompensate)
        CompensateButton.grid(column = 0, row = 3)
        Master.bind('<c>', lambda event:self.SwitchCompensate(True))

    def SwitchCompensate(self, FromBinding = False):
        if FromBinding:
            self.CompensateVar.set(1-self.CompensateVar.get())
        self.Compensate[self.SharedData['CurrentStream']] = bool(self.CompensateVar.get())

    def AddStreamVars(self, Name, Geometry):
        self.StreamsShapes[Name] = Geometry
        self.DisparitiesMaps[Name] = np.zeros(tuple(Geometry[:2]) + (2,))
        self.Compensate[Name] = False

    def Decrypt(self, Socket, t, eventList):
        self.DisparitiesMaps[Socket][eventList[0][0], eventList[0][1], :] = [eventList[1] * eventList[2], t]

    def DelStreamVars(self, Name):
        del self.DisparitiesMaps[Name]
        del self.StreamsShapes[Name]
        del self.Compensate[Name]

    def Update(self, Reload = False):
        if self.Active:
            CurrentStream = self.SharedData['CurrentStream']
            self.CompensateVar.set(self.Compensate[CurrentStream])
            Tau = self.SharedData['Tau']
            DispMap = self.DisparitiesMaps[CurrentStream]
            StreamTime = self.SharedData['t']
            if self.Compensate[CurrentStream]:
                Map = np.zeros(DispMap.shape[:2])
                xs, ys = np.where((StreamTime - DispMap[:,:,1]) < Tau)
                ds = DispMap[:,:,0][xs, ys]
                xs += ds.astype(int)
                Map[xs, ys] = abs(ds)
            else:
                Map = abs(DispMap[:,:,0]) * ((StreamTime - DispMap[:,:,1]) < Tau)

            if Reload:
                self.DisplayImShow = self.Ax.imshow(np.transpose(Map), vmin = self.MinDisparity, vmax = self.MaxDisparity, origin = "lower", cmap = 'hot')
                self.Ax.set_xlim(-0.5, self.StreamsShapes[CurrentStream][0]-0.5)
                self.Ax.set_ylim(-0.5, self.StreamsShapes[CurrentStream][1]-0.5)
            else:
                self.DisplayImShow.set_data(np.transpose(Map))
            self.EventsLabel['text'] = "{0}".format(int(Map.sum()))
            self.Canvas.draw()
        else:
            if Reload:
                self.EventsLabel['text'] = "0"

    def ChangeMaxD(self, var):
        self.MaxDisparity = max(self.MinDisparity + 5, self.MaxDisparity + var * 5)
        self.MaxDisparityLabel['text'] = "Max Disparity : {0}".format(self.MaxDisparity)
        if self.Active:
            self.DisplayImShow.set_clim((self.MinDisparity, self.MaxDisparity))
        self.DrawColorbar()

    def ChangeMinD(self, var):
        self.MinDisparity = max(0, min(self.MaxDisparity-5, self.MinDisparity + var * 5))
        self.MinDisparityLabel['text'] = "Min Disparity : {0}".format(self.MinDisparity)
        if self.Active:
            self.DisplayImShow.set_clim((self.MinDisparity, self.MaxDisparity))
        self.DrawColorbar()

    def AutoSet(self):
        DispMap = self.DisparitiesMaps[self.SharedData['CurrentStream']]
        ds = abs(DispMap[:,:,0][np.where((self.SharedData['t'] - DispMap[:,:,1]) < self.SharedData['Tau'])])
        if ds.shape[0] == 0 or int(ds.max()) == 0:
            self.MinDisparity = 0
            self.MaxDisparity = 50
        else:
            self.MinDisparity = int(ds.min()) - (int(ds.min())%5)
            self.MaxDisparity = int(ds.max()) + 5 - int(ds.max())%5
            N = ds.shape[0]
            while True:
                if (ds < (self.MinDisparity + 5)).sum() < 0.01 * N:
                    self.MinDisparity = self.MinDisparity + 5
                else:
                    break
            while self.MaxDisparity >= self.MinDisparity + 5:
                if (ds > (self.MaxDisparity - 5)).sum() < 0.01 * N:
                    self.MaxDisparity = self.MaxDisparity - 5
                else:
                    break
        self.MinDisparityLabel['text'] = "Min Disparity : {0}".format(self.MinDisparity)
        self.MaxDisparityLabel['text'] = "Max Disparity : {0}".format(self.MaxDisparity)
        if self.Active:
            self.DisplayImShow.set_clim((self.MinDisparity, self.MaxDisparity))
        self.DrawColorbar()

    def AutoSetOld(self):
        DispMap = self.DisparitiesMaps[self.SharedData['CurrentStream']]
        ds = abs(DispMap[:,:,0][np.where((self.SharedData['t'] - DispMap[:,:,1]) < self.SharedData['Tau'])])
        if ds.shape[0] == 0 or int(ds.max()) == 0:
            self.MinDisparity = 0
            self.MaxDisparity = 30
        else:
            self.MinDisparity = int(ds.min())
            self.MaxDisparity = int(ds.max())
        if self.MinDisparity%5 != 0:
            self.MinDisparity = int(self.MinDisparity - self.MinDisparity%5)
        if self.MaxDisparity%5 != 0:
            self.MaxDisparity = int(self.MaxDisparity - self.MaxDisparity%5 + 5)
        self.MinDisparityLabel['text'] = "Min Disparity : {0}".format(self.MinDisparity)
        self.MaxDisparityLabel['text'] = "Max Disparity : {0}".format(self.MaxDisparity)
        if self.Active:
            self.DisplayImShow.set_clim((self.MinDisparity, self.MaxDisparity))
        self.DrawColorbar()

    def DrawColorbar(self):
        self.CBAx.cla()
        norm = matplotlib.colors.Normalize(vmin=self.MinDisparity, vmax=self.MaxDisparity)

        self.CB = matplotlib.colorbar.ColorbarBase(self.CBAx, cmap=plt.cm.hot,
                                norm=norm,
                                orientation='horizontal')
        self.CBCanvas.draw()

class TrackerHandler:
    Key = 2
    DisplayType = "Dot"
    Name = "Trackers"
    def __init__(self, Master, Display, Canvas, InfoFrame, OptionsFrame, SharedData):
        self.Ax = Display
        self.Canvas = Canvas
        self.Active = ((self.Key == 1) or (self.DisplayType == "Dot"))
        self.SharedData = SharedData
        self.Reload = SharedData['ReloadCommand']

        self.TrackersScalingSize = 20
        self.ActiveTrackers = {'None': {}}
        self.UpdatedTrackers = {}
        self.DisplayedTrackers = {}
        self.DisplayedTrackersIDs = {}
        self.CurrentStreamLastUpdatedAlphasTs = {}
        L = Tk.Label(InfoFrame, text = "Active trackers :")
        L.pack(anchor = Tk.W)
        self.TrackersLabel = Tk.Label(InfoFrame)
        self.TrackersLabel.pack(anchor = Tk.W)
        
    def Decrypt(self, Socket, t, eventList):
        ID = eventList.pop(1)
        self.ActiveTrackers[Socket][ID] = eventList + [t,t]
        self.UpdatedTrackers[Socket].add(ID)
    def AddStreamVars(self, Name, Geometry):
        self.ActiveTrackers[Name] = {}
        self.UpdatedTrackers[Name] = set()
        self.CurrentStreamLastUpdatedAlphasTs[Name] = -np.inf
    def DelStreamVars(self, Name):
        del self.ActiveTrackers[Name]
        del self.UpdatedTrackers[Name]
        del self.CurrentStreamLastUpdatedAlphasTs[Name]
    
    def AddTrackerParts(self, TrackerID, TrackerData, Alpha):
        self.DisplayedTrackers[TrackerID] = [None, None]
        self.DisplayedTrackers[TrackerID][0] = self.Ax.plot(TrackerData[0][0], TrackerData[0][1], color = TrackerData[3], marker = TrackerData[4], alpha = Alpha)[0]
        Theta, S = TrackerData[1], TrackerData[2]
        dx, dy = self.TrackersScalingSize * S * np.cos(Theta), self.TrackersScalingSize * S * np.sin(Theta)
        self.DisplayedTrackers[TrackerID][1] = self.Ax.plot([TrackerData[0][0], TrackerData[0][0]+dx], [TrackerData[0][1], TrackerData[0][1]+dy], color = TrackerData[3], alpha = Alpha)[0]
        self.DisplayedTrackersIDs[TrackerID] = self.Ax.text(TrackerData[0][0] + 5, TrackerData[0][1], str(TrackerID), color = TrackerData[3], alpha = Alpha)

    def UpdateTrackerParts(self, TrackerID, TrackerData):
        self.DisplayedTrackers[TrackerID][0].set_data(TrackerData[0][0], TrackerData[0][1])
        self.DisplayedTrackers[TrackerID][0].set_color(TrackerData[3])
        self.DisplayedTrackers[TrackerID][0].set_marker(TrackerData[4])
        Theta, S = TrackerData[1], TrackerData[2]
        dx, dy = self.TrackersScalingSize * S * np.cos(Theta), self.TrackersScalingSize * S * np.sin(Theta)
        self.DisplayedTrackers[TrackerID][1].set_data([TrackerData[0][0], TrackerData[0][0]+dx], [TrackerData[0][1], TrackerData[0][1]+dy])
        self.DisplayedTrackers[TrackerID][1].set_color(TrackerData[3])
        self.DisplayedTrackersIDs[TrackerID].set_x(TrackerData[0][0] + 5)
        self.DisplayedTrackersIDs[TrackerID].set_y(TrackerData[0][1])
        self.DisplayedTrackersIDs[TrackerID].set_color(TrackerData[3])

    def UpdateTrackerAlpha(self, TrackerID, Alpha):
        self.DisplayedTrackers[TrackerID][0].set_alpha(Alpha)
        self.DisplayedTrackers[TrackerID][1].set_alpha(Alpha)
        self.DisplayedTrackersIDs[TrackerID].set_alpha(Alpha)
    def RemoveTrackerParts(self, TrackerID):
        for Part in self.DisplayedTrackers[TrackerID]:
            Part.remove()
        self.DisplayedTrackersIDs[TrackerID].remove()
        del self.DisplayedTrackers[TrackerID]
        del self.DisplayedTrackersIDs[TrackerID]

    def Update(self, Reload = False):
        CurrentStream = self.SharedData['CurrentStream']
        Tau = self.SharedData['Tau']
        if self.Active:
            T = self.SharedData['t']
            if Reload:
                self.DisplayedTrackers = {}
                self.DisplayedTrackersIDs = {}
                for TrackerID, TrackerData in dict(self.ActiveTrackers[CurrentStream]).items():
                    Alpha = min(1, max(0, np.e**((TrackerData[-2] - T)/Tau)))
                    self.AddTrackerParts(TrackerID, TrackerData, Alpha)
                self.UpdatedTrackers[CurrentStream].clear()
            else:
                UpdateAlphas = (T - self.CurrentStreamLastUpdatedAlphasTs[CurrentStream] > Tau / 5)
                if UpdateAlphas:
                    self.CurrentStreamLastUpdatedAlphasTs[CurrentStream] = T
                    for TrackerID in list(self.DisplayedTrackers.keys()):
                        TrackerData = self.ActiveTrackers[CurrentStream][TrackerID]
                        Alpha = min(1, max(0, np.e**((TrackerData[-2] - T)/Tau)))
                        if Alpha < 0.001:
                            self.RemoveTrackerParts(TrackerID)
                            continue
                        if TrackerID in self.UpdatedTrackers[CurrentStream]:
                            if TrackerID not in self.DisplayedTrackers.keys():
                                self.AddTrackerParts(TrackerID, TrackerData, Alpha)
                            else:
                                self.UpdateTrackerParts(TrackerID, TrackerData)
                            self.UpdateTrackerAlpha(TrackerID, Alpha)
                else:
                    for TrackerID in list(self.UpdatedTrackers[CurrentStream]):
                        TrackerData = self.ActiveTrackers[CurrentStream][TrackerID]
                        Alpha = min(1, max(0, np.e**((TrackerData[-2] - T)/Tau)))
                        if TrackerID not in self.DisplayedTrackers.keys():
                            self.AddTrackerParts(TrackerID, TrackerData, Alpha)
                        else:
                            if Alpha < 0.001:
                                self.RemoveTrackerParts(TrackerID)
                            else:
                                self.UpdateTrackerParts(TrackerID, TrackerData)
                                self.UpdateTrackerAlpha(TrackerID, Alpha)
                self.UpdatedTrackers[CurrentStream].clear()
        self.TrackersLabel['text'] = "{0}".format(len(self.ActiveTrackers[CurrentStream]))

class FlowHandler:
    Key = 6
    DisplayType = "Dot"
    Name = "Optical Flow"
    def __init__(self, Master, Display, Canvas, InfoFrame, OptionsFrame, SharedData):
        self.Ax = Display
        self.Canvas = Canvas
        self.Active = ((self.Key == 1) or (self.DisplayType == "Dot"))
        self.SharedData = SharedData
        self.Reload = SharedData['ReloadCommand']

        self.AlphaPerPixel = 0.4
        self.PlottedFlowDensity = 0.1
        self.LastClearAndUpdate = -np.inf
        self.FlowMaps = {}
        self.UpdatedFlowsLocations = {}
        self.PlottedFlows = {}
        L = Tk.Label(InfoFrame, text = "Plotted Flows :")
        L.pack(anchor = Tk.W)
        self.FlowsLabel = Tk.Label(InfoFrame)
        self.FlowsLabel.pack(anchor = Tk.W)
        
        self.FDLabel = Tk.Label(OptionsFrame, text = "Flow density : {0:.1f}".format(self.PlottedFlowDensity))
        self.FDLabel.grid(column = 0, row = 0)
        FDMinusButton = Tk.Button(OptionsFrame, text = ' - ', command = lambda: self.ChangeFD(-1))
        FDMinusButton.grid(column = 1, row = 0)
        FDPlusButton = Tk.Button(OptionsFrame, text = ' + ', command = lambda: self.ChangeFD(+1))
        FDPlusButton.grid(column = 2, row = 0)

    def ChangeFD(self, var):
        self.PlottedFlowDensity = max(0, min(1., self.PlottedFlowDensity + var * 0.1))
        self.FDLabel['text'] = "Flow density : {0:.1f}".format(self.PlottedFlowDensity)

    def Decrypt(self, Socket, t, eventList):
        location, flow = eventList
        if (abs(flow[0]) - int(abs(flow[0]))) >= self.PlottedFlowDensity:
            return
        self.FlowMaps[Socket][location[0], location[1], :2] = np.array(flow)
        self.FlowMaps[Socket][location[0], location[1], 2] = t
        self.UpdatedFlowsLocations[Socket].add(tuple(location))
    def AddStreamVars(self, Name, Geometry):
        self.FlowMaps[Name] = np.zeros(tuple(Geometry[:2]) + (3,))
        self.FlowMaps[Name][:,:,2] = -np.inf
        self.UpdatedFlowsLocations[Name] = set()
        self.PlottedFlows[Name] = []

    def DelStreamVars(self, Name):
        del self.FlowMaps[Name]
        del self.UpdatedFlowsLocations[Name]
        del self.PlottedFlows[Name]
    
    def PlotFlow(self, x, y, T, Tau, CurrentStream):
        fx, fy, t = self.FlowMaps[CurrentStream][x, y, :]
        f = np.array([fx, fy])
        n = np.linalg.norm(f)
        alpha = self.alpha(t, T, Tau)
        if alpha > 0.1:
            nf = f / n
            r = max(0., nf[0])
            g = max(0., (nf * np.array([-1, np.sqrt(3)])).sum()/2)
            b = max(0., (nf * np.array([-1, -np.sqrt(3)])).sum()/2)
            color = np.sqrt(np.minimum(np.array([r, g, b, alpha]), np.ones(4, dtype = float)))
            df = f * Tau
            self.PlottedFlows[CurrentStream] += [(t, n, self.Ax.arrow(x, y, df[0], df[1], color = color))]

    def alpha(self, t, T, Tau):
        return max(0., 1-(T-t)/Tau)

    def Update(self, Reload = False):
        CurrentStream = self.SharedData['CurrentStream']
        Tau = self.SharedData['Tau']
        if self.Active:
            T = self.SharedData['t']
            if Reload:
                self.LastClearAndUpdate = T
                for t, n, arrow in self.PlottedFlows[CurrentStream]:
                    arrow.remove()
                self.PlottedFlows[CurrentStream] = []
                for x in range(self.FlowMaps[CurrentStream].shape[0]):
                    for y in range(self.FlowMaps[CurrentStream].shape[1]):
                        self.PlotFlow(x, y, T, Tau, CurrentStream)
            else:
                if T - self.LastClearAndUpdate > Tau / 5:
                    KeptArrows = []
                    for n, (t, n, arrow) in enumerate(self.PlottedFlows[CurrentStream]):
                        alpha = self.alpha(t, T, Tau)
                        if alpha < 0.1:
                            arrow.remove()
                        else:
                            arrow.set_alpha(alpha)
                            KeptArrows += [(t, n, arrow)]
                    self.PlottedFlows[CurrentStream] = KeptArrows
                    LastClearAndUpdate = T
                for x, y in set(self.UpdatedFlowsLocations[CurrentStream]):
                    self.PlotFlow(x, y, T, Tau, CurrentStream)
                self.UpdatedFlowsLocations[CurrentStream].clear()

        self.FlowsLabel['text'] = "{0}".format(len(self.PlottedFlows[CurrentStream]))


_HANDLERS = [EventHandler, TrackerHandler, DisparityHandler, FlowHandler]

class Display:
    _DefaultTau = 0.030
    _AutoRecordDt = 0.01
    _AutoRecordFolder = 'home/dardelet/Pictures/Recordings/'
    def __init__(self, mainPort = 54242, questionPort = 54243, responsePort = 54244, listen = "", responseAddress = 'localhost'):
        self.SharedData = {'t':0, 'Tau':self._DefaultTau, 'Stream':None, 'ReloadCommand': self.Reload}
        self.dTau = 0.005
        self.PreviousTime = 0.
        self.PreviousNEvents = 0
	
        self.MainSocket = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
        self.MainSocket.bind((listen, mainPort))

        self.QuestionSocket = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
        self.QuestionSocket.bind((listen, questionPort))

        self.ResponseSocket = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
        ResponseAddress = ("localhost", responsePort)

        self.CurrentlyDisplayedStreams = []
        self.Streams = []
        self.StreamsTypes = {}
        self.LastTSAtRTDUpdate = {}
        self.StreamsInfos = {}
        self.StreamsTaus = {}
        self.StreamsTimes = {}

        self.ReloadNames = False
        self.ReloadMap = True
        
        self.MainWindow = Tk.Tk()
        self.MainWindow.title('Projector')
        # Create different frames
        self.DisplayFrame = Tk.Frame(self.MainWindow)
        self.DisplayFrame.grid(row = 0, column = 0, rowspan = 8)

        self.RTFrame = Tk.Frame(self.MainWindow)
        self.RTFrame.grid(row = 0, column = 2, rowspan = 2, sticky = Tk.NW)

        self.HandlersInfoFrame = Tk.Frame(self.MainWindow)
        self.HandlersInfoFrame.grid(row = 1, column = 1, columnspan = 2, sticky = Tk.W)
        self.HandlersMapInfoFrame = Tk.Frame(self.HandlersInfoFrame)
        self.HandlersMapInfoFrame.pack(anchor = Tk.W)
        self.HandlersDotInfoFrame = Tk.Frame(self.HandlersInfoFrame)
        self.HandlersDotInfoFrame.pack(anchor = Tk.W)

        self.DisplayedTypesFrame = Tk.Frame(self.MainWindow)
        self.DisplayedTypesFrame.grid(row = 3, column = 1, columnspan = 2, sticky = Tk.EW)

        self.StreamsFrame = Tk.Frame(self.MainWindow)
        self.StreamsFrame.grid(row = 4, column = 1, columnspan = 2, sticky = Tk.W)

        self.StreamInfoFrame = Tk.Frame(self.MainWindow)
        self.StreamInfoFrame.grid(row = 0, column = 1, sticky = Tk.W)

        TauFrame = Tk.Frame(self.StreamInfoFrame)
        TauFrame.pack(anchor = Tk.W)

        self.Display = Figure(figsize=(7,7), dpi=150)
        self.DisplayAx = self.Display.add_subplot(111)
        DisplayCanvas = FigureCanvasTkAgg(self.Display, self.DisplayFrame)
        DisplayCanvas.mpl_connect('button_press_event', self.OnClick)
        DisplayCanvas.draw()
        DisplayCanvas.get_tk_widget().grid(row = 0, column = 0)
        

        L = Tk.Label(self.DisplayedTypesFrame, text = "Display selection :")
        L.grid(row = 0, column = 0, sticky = Tk.W)

        self.DisplayedTypesMapFrame = Tk.Frame(self.DisplayedTypesFrame, borderwidth = 1, relief=Tk.SOLID)
        self.DisplayedTypesMapFrame.grid(row = 1, column = 0, sticky = Tk.EW)
        self.DisplayedTypesDotFrame = Tk.Frame(self.DisplayedTypesFrame,  borderwidth = 1, relief=Tk.SOLID)
        self.DisplayedTypesDotFrame.grid(row = 2, column = 0, sticky = Tk.EW)
        self.DisplayedTypesVars = {}
        self.DisplayedMapTypeVar = Tk.IntVar(master = self.MainWindow)
        self.DisplayedMapTypeVar.set(1)
        self.Handlers = {}
        UsedKeys = set()
        HandlersRows = {"Map":0, "Dot":0}
        for nType, Handler in enumerate(_HANDLERS):
            letter = Handler.Name[0].lower()
            while letter in UsedKeys:
                letter = alphabet[alphabet.index(letter)-1]
            UsedKeys.add(letter)
            if Handler.DisplayType == "Dot":
                HandlerFrame = self.DisplayedTypesDotFrame
                HandlerInfoFrame = self.HandlersDotInfoFrame
                self.DisplayedTypesVars[Handler.Key] = Tk.IntVar(master = self.MainWindow)
                self.DisplayedTypesVars[Handler.Key].set(1)
                Button = Tk.Checkbutton(HandlerFrame, text = Handler.Name + ' ({0})'.format(letter), variable = self.DisplayedTypesVars[Handler.Key], command = lambda E=Handler.Key : self.SwitchDisplayedTypeOnOff(E))
                self.MainWindow.bind('<{0}>'.format(letter), lambda event, E = Handler.Key:self.SwitchDisplayedTypeOnOff(E, True))
            elif Handler.DisplayType == "Map":
                HandlerFrame = self.DisplayedTypesMapFrame
                HandlerInfoFrame = self.HandlersMapInfoFrame
                Button = Tk.Radiobutton(HandlerFrame, text = Handler.Name + ' ({0})'.format(letter), variable = self.DisplayedMapTypeVar, value = Handler.Key, command = lambda:self.SwitchDisplayedMapType())
                self.MainWindow.bind('<{0}>'.format(letter), lambda event, E = Handler.Key:self.SwitchDisplayedMapType(E))
            else:
                raise Exception("Ill defined handler")
            Button.grid(row = HandlersRows[Handler.DisplayType], column = 0, sticky = Tk.W)
            HandlerOptionsFrame = Tk.Frame(HandlerFrame)
            HandlerOptionsFrame.grid(row = HandlersRows[Handler.DisplayType], column = 1, sticky = Tk.W)
            HandlersRows[Handler.DisplayType] += 1

            HandlerInfoFrame = Tk.Frame(HandlerInfoFrame)
            HandlerInfoFrame.pack(anchor = Tk.W)
            
            self.Handlers[Handler.Key] = Handler(self.MainWindow, self.DisplayAx, self.Display.canvas, HandlerInfoFrame, HandlerOptionsFrame, self.SharedData)
        for UnusedKey in range(10):
            if UnusedKey not in self.Handlers:
                self.Handlers[UnusedKey] = NullHandler()

        self.StreamsVariable = Tk.IntVar(master = self.MainWindow)
        self._MinimumGeometry = np.array([10, 10, 2])
        self.AddStream('None', self._MinimumGeometry)

        self.RTMinLim = 10**-4
        self.RTMaxLim = 10**1
        self.CurrentRTValue = 1
        self.UpdateRTDTimeConstant = 500 # in ms
        self.RealTimeDisplay = Figure(figsize=(1,2.5), dpi=100)
        self.RealTimeDisplayAx = self.RealTimeDisplay.add_subplot(111)
        self.RealTimeDisplayAx.tick_params('both', bottom =  'off', labelbottom = 'off')
        self.RealTimeDisplayAx.set_yscale('log')
        LogValue = np.log10(self.CurrentRTValue)
        self.RealTimeDisplayAx.set_title("{0:.1f}".format(LogValue))
        self.RealTimeDisplayAx.set_ylim(self.RTMinLim, self.RTMaxLim)
        self.RealTimeDisplayAx.set_xlim(0,1)

        self.RTDColors = [C.get_rgb() for C in Color("red").range_to(Color("green"),10)]
        self.RTRectangle = Rectangle((0,0), 1, self.CurrentRTValue)
        self.RTRectangle.set_color(self.RTDColors[max(min(int(10 * (np.log10(self.CurrentRTValue) + 4) / 5), 9), 0)])
        self.RealTimeDisplayAx.add_patch(self.RTRectangle)

        RTDisplayCanvas = FigureCanvasTkAgg(self.RealTimeDisplay, self.RTFrame)
        RTDisplayCanvas.draw()
        self.RealTimeDisplay.tight_layout()
        RTDisplayCanvas.get_tk_widget().grid(row = 0, column = 0, sticky = Tk.NW)

        self.FPSValue = 0.
        self.FPSLabel = Tk.Label(self.RTFrame, text = "FPS : {0}".format(0))
        self.FPSLabel.grid(row = 1, column = 0, sticky = Tk.NW)

        L = Tk.Label(self.StreamsFrame, text = "Streams selection :")
        L.pack(anchor = Tk.W)
        self.StreamsButtons = {'None': Tk.Radiobutton(self.StreamsFrame, text = self.StreamsInfos['None'], variable = self.StreamsVariable, value = 0, command = lambda:self.SwitchStream(0))}
        self.StreamsButtons['None'].pack(anchor = Tk.W)
        self.SharedData['CurrentStream'] = self.Streams[self.StreamsVariable.get()]

        L = Tk.Label(self.StreamInfoFrame, text = "Current Timestamp :")
        L.pack(anchor = Tk.W)
        self.TSLabel = Tk.Label(self.StreamInfoFrame)
        self.TSLabel.pack(anchor = Tk.W)
        L = Tk.Label(self.StreamInfoFrame, text = "Pipeline efficiency :")
        L.pack(anchor = Tk.W)
        self.EfficiencyLabel = Tk.Label(self.StreamInfoFrame)
        self.EfficiencyLabel.pack(anchor = Tk.W)
        L = Tk.Label(self.StreamInfoFrame, text = "Point location :")
        L.pack(anchor = Tk.W)
        self.PointLabel = Tk.Label(self.StreamInfoFrame)
        self.PointLabel.pack(anchor = Tk.W)

        self.TauLabel = Tk.Label(TauFrame, text = "Tau : {0} ms".format(int(1000*self.Tau)))
        self.TauLabel.grid(column = 0, row = 0)
        TauMinusButton = Tk.Button(TauFrame, text = ' - ', command = lambda: self.ChangeTau(-1))
        TauMinusButton.grid(column = 1, row = 0)
        TauPlusButton = Tk.Button(TauFrame, text = ' + ', command = lambda: self.ChangeTau(+1))
        TauPlusButton.grid(column = 2, row = 0)
        self.dTauLabel = Tk.Label(TauFrame, text = "dTau : {0:.2f} ms".format(int(1000*self.dTau)))
        self.dTauLabel.grid(column = 0, row = 1)
        dTauDivideButton = Tk.Button(TauFrame, text = ' /2 ', command = lambda: self.ChangedTau(-1))
        dTauDivideButton.grid(column = 1, row = 1)
        dTauMultiplyButton = Tk.Button(TauFrame, text = ' *2 ', command = lambda: self.ChangedTau(+1))
        dTauMultiplyButton.grid(column = 2, row = 1)

        self.MainWindow.bind('<Escape>', lambda event: self._on_closing())
        self.MainWindow.bind('<KP_Add>', lambda event: self.ChangeTau(+1))
        self.MainWindow.bind('<KP_Subtract>', lambda event: self.ChangeTau(-1))
        self.MainWindow.bind('<KP_Multiply>', lambda event: self.ChangedTau(+1))
        self.MainWindow.bind('<KP_Divide>', lambda event: self.ChangedTau(-1))
        self.MainWindow.bind('<KP_0>', lambda event: self.SwitchDisplayedPolasOnOff(0,0))
        self.MainWindow.bind('<KP_1>', lambda event: self.SwitchDisplayedPolasOnOff(1,1))
        self.MainWindow.bind('<KP_2>', lambda event: self.SetPolasMode(2))
        self.MainWindow.bind('<Down>', lambda event: self.SwitchStream(+1))
        self.MainWindow.bind('<Up>', lambda event: self.SwitchStream(-1))
        self.MainWindow.protocol('WM_DELETE_WINDOW', self._on_closing)

        ClearStreamsButton = Tk.Button(self.MainWindow, text='Clear', command = self.ClearStreams)
        ClearStreamsButton.grid(row = 5, column = 1, sticky = Tk.W)


        self.QuestionThread = QuestionThreadClass(self, self.QuestionSocket)
        self.QuestionThread.start()

        self.ResponseBot = ResponseBotClass(self, self.ResponseSocket, ResponseAddress)

        self.MainThread = MainThreadClass(self, self.MainSocket)
        self.MainThread.start()

        self.UpdateStreamsList()
        self.UpdateScene()
        self.UpdateTS()
        self.UpdateRTD()
        self.MainWindow.mainloop()

    @property
    def Tau(self):
        return self.SharedData['Tau']
    @property
    def CurrentStream(self):
        return self.SharedData['CurrentStream']

    def AddStream(self, Name, Geometry, StreamInfo = "No display", StreamType = 'Stream'):
        Geometry = np.maximum(Geometry, self._MinimumGeometry)
        self.Streams += [Name]
        self.StreamsTypes[Name] = StreamType
        self.StreamsTimes[Name] = 0
        self.StreamsTaus[Name] = self._DefaultTau
        self.LastTSAtRTDUpdate[Name] = 0
        self.StreamsInfos[Name] = StreamInfo

        for Handler in self.Handlers.values():
            Handler.AddStreamVars(Name, Geometry)

    def RemoveStream(self, Name):
        del self.StreamsTimes[Name]
        del self.LastTSAtRTDUpdate[Name]
        del self.StreamsInfos[Name]
        del self.StreamsTypes[Name]
        del self.StreamsTaus[Name]
        for Handler in self.Handlers.values():
            Handler.DelStreamVars(Name)
        self.CurrentlyDisplayedStreams.remove(Name)
        self.Streams.remove(Name)
        self.SharedData['CurrentStream'] = self.Streams[-1]

    def OnClick(self, event):
        if event.inaxes is None:
            return
        self.PointLabel['text'] = "x = {0}, y = {1}".format(int(event.xdata), int(event.ydata))

    def SwitchStream(self, var):
        if var != 0:
            self.StreamsVariable.set(min(max(self.StreamsVariable.get() + var, 0), len(self.Streams)-1))
        self.Reload()

    def SwitchDisplayedMapType(self, NewValue = None):
        if NewValue is None:
            NewValue = self.DisplayedMapTypeVar.get()
        else:
            self.DisplayedMapTypeVar.set(NewValue)
        for Handler in self.Handlers.values():
            if Handler.DisplayType == 'None':
                continue
            if Handler.DisplayType == "Map":
                Handler.Active = (Handler.Key == NewValue)
        self.Reload()

    def SwitchDisplayedTypeOnOff(self, EventKey, FromBinding = False):
        if FromBinding:
            self.DisplayedTypesVars[EventKey].set(1 - self.DisplayedTypesVars[EventKey].get())
        self.Handlers[EventKey].Active = bool(self.DisplayedTypesVars[EventKey].get())
        self.Reload()

    def ClearStreams(self):
        Except = ['None']
        if self.StreamsVariable.get() > 0:
            Except += [self.Streams[self.StreamsVariable.get()]]
        PreviousStreams = list(self.Streams)
        for Stream in PreviousStreams:
            if not Stream in Except:
                self.RemoveStream(Stream)
        self.ReloadNames = True
        self.UpdateStreamsList()

    def UpdateStreamsList(self):
        if self.CurrentlyDisplayedStreams != self.Streams or self.ReloadNames:
            self.ReloadNames = False
            self.CurrentlyDisplayedStreams = []
            print("Updating stream list, aiming for {0} streams".format(len(self.Streams)))
            for Stream in dict(self.StreamsButtons).keys():
                self.StreamsButtons[Stream].destroy()
                del self.StreamsButtons[Stream]
            for Stream in self.Streams:
                self.StreamsButtons[Stream] = Tk.Radiobutton(self.StreamsFrame, text = self.StreamsInfos[Stream], variable = self.StreamsVariable, value = self.Streams.index(Stream), command = self.Reload)
                self.StreamsButtons[Stream].pack(anchor = Tk.W)
                self.CurrentlyDisplayedStreams += [Stream]
            self.StreamsVariable.set(len(self.Streams)-1)
            self.Reload()
        self.MainWindow.after(100, self.UpdateStreamsList)

    def ChangeTau(self, var):
        self.StreamsTaus[self.CurrentStream] = max(self.dTau, self.Tau + var * self.dTau)
        self.UpdateTauDisplayed()

    def ChangedTau(self, var):
        self.dTau *= 2**var
        DisplayeddTau = self.dTau
        Prefixes = ['', 'm', 'µ', 'n']
        nPrefix = 0
        while DisplayeddTau < 1 and nPrefix != 3:
            DisplayeddTau *= 1000
            nPrefix += 1
        self.dTauLabel['text'] = "dTau : {0:.1f} {1}s".format(DisplayeddTau, Prefixes[nPrefix])

    def UpdateTauDisplayed(self):
        self.SharedData['Tau'] = self.StreamsTaus[self.CurrentStream]
        DisplayedTau = self.Tau
        Prefixes = ['', 'm', 'µ', 'n']
        nPrefix = 0
        while DisplayedTau < 1 and nPrefix != 3:
            DisplayedTau *= 1000
            nPrefix += 1
        self.TauLabel['text'] = "Tau : {0:.1f} {1}s".format(DisplayedTau, Prefixes[nPrefix])

    def printPolaritiesVariable(self):
        print(self.PolaritiesVariable.get())

    def Reload(self):
        self.SharedData['CurrentStream'] = self.Streams[self.StreamsVariable.get()]
        self.ReloadMap = True

    def _on_closing(self):
        KillDisplay()
        self.MainSocket.close()
        self.QuestionSocket.close()
        self.ResponseSocket.close()

        self.MainWindow.quit()
        self.MainWindow.destroy()

    def UpdateRTD(self):
        TimeEllapsed = self.StreamsTimes[self.CurrentStream] - self.LastTSAtRTDUpdate[self.CurrentStream]
        self.LastTSAtRTDUpdate[self.CurrentStream] = self.StreamsTimes[self.CurrentStream]
        self.CurrentRTValue = max(self.RTMinLim, TimeEllapsed / self.UpdateRTDTimeConstant)

        LogValue = np.log10(self.CurrentRTValue)
        self.RTRectangle.set_color(self.RTDColors[max(min(int(10 * (LogValue + 4) / 5), 9), 0)])
        self.RTRectangle.set_height(self.CurrentRTValue)
        self.RealTimeDisplayAx.set_title("{0:.1f}".format(LogValue))

        self.RealTimeDisplay.canvas.draw()
        self.MainWindow.after(self.UpdateRTDTimeConstant, self.UpdateRTD)

    def UpdateScene(self):
        if self.ReloadMap:
            self.DisplayAx.clear()
        self.SharedData['t'] = self.StreamsTimes[self.CurrentStream]
        for Handler in self.Handlers.values():
            Handler.Update(self.ReloadMap)

        self.Display.canvas.draw()
        self.ReloadMap = False

        t = time.time()
        self.EfficiencyLabel['text'] = "{0} ev/s".format(int(self.PreviousNEvents / (t - self.PreviousTime)))
        self.FPSValue = self.FPSValue * 0.8 + 0.2/(t - self.PreviousTime)
        self.FPSLabel['text'] = "FPS : {0}".format(int(self.FPSValue))
        self.PreviousTime = t
        self.PreviousNEvents = 0
        self.MainWindow.after(5, self.UpdateScene)

    def _DecryptFullEvent(self, Socket, eventList):
        t = eventList[0]
        self.StreamsTimes[Socket] = max(self.StreamsTimes[Socket], t)
        for Extension in eventList[1:]:
            self.Handlers[Extension[0]].Decrypt(Socket, t, Extension[1:])

    def _InitTauEventsVars(self):
        self._EventDecryptFunctions[5] = self._DecryptTauEvent

    def _DecryptTauEvent(self, Socket, t, eventList):
        self.StreamsTaus[Socket] = eventList[0]
        self.UpdateTauDisplayed()
        
    def UpdateTS(self):
        self.TSLabel['text'] = "{0} ms".format(int(1000*self.StreamsTimes[self.CurrentStream]))
        self.MainWindow.after(50, self.UpdateTS)

class QuestionThreadClass(threading.Thread):
    def __init__(self, parent, conn):
        self.Connexion = conn
        self.Parent = parent
        threading.Thread.__init__(self)
        self.Run = True

    def run(self):
        while self.Run:
            data_raw, addr = self.Connexion.recvfrom(1024)
            
            data = cPickle.loads(data_raw)
            print("Received {0} command".format(data['command']))
            print(data)
            if data['command'] == 'isup':
                self.Parent.ResponseBot.Answer(data['id'],True)
            if data['command'] == 'kill':
                self.Run = False
            if data['command'] == 'asksocket':
                idgiven = str(addr[1])
                if idgiven not in self.Parent.Streams:
                    self.Parent.ResponseBot.Answer(data['id'], idgiven)
                    self.Parent.AddStream(idgiven, data['shape'], "Unknown ({0})".format(idgiven))
                else:
                    self.Parent.ResponseBot.Answer(data['id'], 'socketexists')

            if data['command'] == 'askspecificsocket':
                socketasked = data['socket']
                self.Parent.ResponseBot.Answer(data['id'], socketasked)
                if socketasked not in self.Parent.Streams:
                    self.Parent.AddStream(socketasked, data['shape'], "Unknown ({0})".format(socketasked))

            if data['command'] == "cleansocket":
                s = data['socket']
                if s in self.Parent.StreamsTimes.keys():
                    self.Parent.StreamsTimes[s] = -0
                    self.Parent.LastTSAtRTDUpdate[s] = -0
                    self.Parent.ResponseBot.Answer(data['id'],'socketcleansed')
                else:
                    self.Parent.ResponseBot.Answer(data['id'],'cleaningfailed')

            if data['command'] == 'socketdata':
                s = data['socket']
                if s in self.Parent.Streams:
                    self.Parent.StreamsInfos[s] = data['infosline1']+"\n"+data['infosline2']
                    self.Parent.ResponseBot.Answer(data['id'],'datareceived')
                    self.Parent.ReloadNames = True
                else:
                    self.Parent.ResponseBot.Answer(data['id'],'datatransmissionfailed')

            if data['command'] == 'destroysocket':
                s = data['socket']
                try:
                    self.Parent.RemoveStream(s)
                    self.Parent.ResponseBot.Answer(data['id'],'socketdestroyed')
                except:
                    self.Parent.ResponseBot.Answer(data['id'],'destructionfailed')

            if data['command'] == 'rewind':
                Stream = data['socket']
                tNew = data['tNew']
                self.Parent.StreamsTimes[Stream] = tNew
                self.Parent.ResponseBot.Answer(data['id'],'rewinded')

        self.Connexion.close()
        print("Questions connexion closed")

class ResponseBotClass:
    def __init__(self, parent, conn, address):
        self.Address = address
        self.Connexion = conn
    
    def Answer(self, questionId, answer):
        answerdict = {'id': questionId, 'answer':answer}
        self.Connexion.sendto(cPickle.dumps(answerdict), self.Address)

class MainThreadClass(threading.Thread):
    def __init__(self, parent, conn):
        self.Connexion = conn
        self.Parent = parent
        threading.Thread.__init__(self)
        self.Run = True
        self.WritingBuffer = {}

    def run(self):
        while self.Run:
            data, addr = self.Connexion.recvfrom(8192)
            if b'kill' in data:
                self.Run = False
            else:
                Update = cPickle.loads(data)
                SocketConcerned = Update[0]

                if SocketConcerned in self.Parent.Streams:
                    for Event in Update[1:]:
                        self.Parent.PreviousNEvents += 1
                        self.Parent._DecryptFullEvent(SocketConcerned, Event)

        self.Connexion.close()
        print("Main connexion closed")

D = Display()

