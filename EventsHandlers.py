import numpy as np
import tkinter as Tk
import random

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class NullHandler:
    DisplayType = 'None'
    def __init__(self):
        pass
    def Decrypt(self, Socket, t, Data, Location):
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

        self.EventsLabel = Tk.Label(InfoFrame, text = "Events shown : 0")
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
        
    def Decrypt(self, Socket, t, Data, Location = None):
        self.StreamsMaps[Socket][Data[0][0], Data[0][1], Data[1]] = t
    def AddStreamVars(self, Name, Geometry):
        self.StreamsShapes[Name] = Geometry
        self.StreamsMaps[Name] = -10*np.ones(tuple(self.StreamsShapes[Name]) + (2,))
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
            self.EventsLabel['text'] = "Events shown : {0}".format(int(Map.sum()))
        else:
            if Reload:
                self.EventsLabel['text'] = "Events shown : 0"

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

        self.DisparitiesLabel = Tk.Label(InfoFrame, text = "Disparities shown : 0")
        self.DisparitiesLabel.grid(row = 0, column = 0, sticky = Tk.NW)

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
        self.DisparitiesMaps[Name] = np.zeros(tuple(Geometry) + (2,))
        self.Compensate[Name] = False

    def Decrypt(self, Socket, t, Data, Location = None):
        if Location is None:
            return
        self.DisparitiesMaps[Socket][Location[0], Location[1], :] = [Data[0] * Data[1], t]

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
            self.DisparitiesLabel['text'] = "Disparities shown : {0}".format(int((Map>0).sum()))
            self.Canvas.draw()
        else:
            if Reload:
                self.DisparitiesLabel['text'] = "Disparities shown : 0"

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
        self.TrackersLabel = Tk.Label(InfoFrame, text = "Active trackers : 0")
        self.TrackersLabel.pack(anchor = Tk.W)
        
    def Decrypt(self, Socket, t, Data, Location = None):
        ID = Data.pop(1)
        self.ActiveTrackers[Socket][ID] = Data + [t,t]
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
        self.TrackersLabel['text'] = "Active trackers : {0}".format(len(self.ActiveTrackers[CurrentStream]))

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
        self.FlowsLabel = Tk.Label(InfoFrame, text = "Plotted Flows : 0")
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

    def Decrypt(self, Socket, t, Data, Location = None):
        if Location is None:
            return
        if random.random() >= self.PlottedFlowDensity:
            return
        self.FlowMaps[Socket][Location[0], Location[1], :2] = np.array(Data)
        self.FlowMaps[Socket][Location[0], Location[1], 2] = t
        self.UpdatedFlowsLocations[Socket].add(tuple(Location))
    def AddStreamVars(self, Name, Geometry):
        self.FlowMaps[Name] = np.zeros(tuple(Geometry) + (3,))
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
            self.FlowsLabel['text'] = "Plotted Flows : {0}".format(len(self.PlottedFlows[CurrentStream]))
        else:
            self.FlowsLabel['text'] = "Plotted Flows : 0"

class TwistHandler:
    Key = 7
    DisplayType = "Info"
    Name = "Twist"
    def __init__(self, Master, Display, Canvas, InfoFrame, OptionsFrame, SharedData):
        self.SharedData = SharedData

        self.ReceivedTwists = {}

        L = Tk.Label(InfoFrame, text = "Latest twist :")
        L.pack(anchor = Tk.W)
        self.TwistLabel = Tk.Label(InfoFrame)
        self.TwistLabel.pack(anchor = Tk.W)        

    def Decrypt(self, Socket, t, Data, Location = None):
        self.ReceivedTwists[Socket] = (Data[0], Data[1])

    def AddStreamVars(self, Name, Geometry):
        self.ReceivedTwists[Name] = (np.zeros(3), np.zeros(3))

    def DelStreamVars(self, Name):
        del self.ReceivedTwists[Name]
    
    def Update(self, Reload = False):
        CurrentStream = self.SharedData['CurrentStream']
        omega, v = self.ReceivedTwists[CurrentStream]
        self.TwistLabel['text'] = "w_x = {0:.3f}rad/s\nw_y = {1:.3f}rad/s\nw_z = {2:.3f}rad/s\nv_x = {3:.3f}m/s\nv_y = {4:.3f}m/s\nv_z = {5:.3f}m/s".format(omega[0], omega[1], omega[2], v[0], v[1], v[2])

_HANDLERS = [EventHandler, TrackerHandler, DisparityHandler, FlowHandler, TwistHandler]

