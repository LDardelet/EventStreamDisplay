import numpy as np
import socket
import matplotlib
import matplotlib.pyplot as pyl
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

import threading
import datetime

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

class Display:
    def __init__(self, mainPort = 54242, questionPort = 54243, responsePort = 54244, listen = "", responseAddress = 'localhost'):
        self.Tau = 30 #in ms
        self.dTau = 5
        self.PreviousTime = 0.
        self.PreviousNEvents = 0
	
        self.MainSocket = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
        self.MainSocket.bind((listen, mainPort))

        self.QuestionSocket = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
        self.QuestionSocket.bind((listen, questionPort))

        self.ResponseSocket = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
        ResponseAddress = ("localhost", responsePort)

        self.QuestionThread = QuestionThreadClass(self, self.QuestionSocket)
        self.QuestionThread.start()

        self.ResponseBot = ResponseBotClass(self, self.ResponseSocket, ResponseAddress)

        self.CurrentlyDisplayedStreams = []
        self.Streams = []
        self.StreamsTypes = {}
        self.LastTSAtRTDUpdate = {}
        self.StreamsInfos = {}
        self.StreamsTimes = {}

        self.ReloadNames = False
        self.ReloadMap = False
        
        self.MainWindow = Tk.Tk()
        self.MainWindow.title('Projector')
        # Create different frames
        self.DisplayFrame = Tk.Frame(self.MainWindow)
        self.DisplayFrame.grid(row = 0, column = 0, rowspan = 8)

        self.RTFrame = Tk.Frame(self.MainWindow)
        self.RTFrame.grid(row = 0, column = 2, rowspan = 2)

        self.EventSpecificFrame = Tk.Frame(self.MainWindow)
        self.EventSpecificFrame.grid(row = 1, column = 1, sticky = Tk.W)

        self.PolaritiesFrame = Tk.Frame(self.MainWindow)
        self.PolaritiesFrame.grid(row = 2, column = 1, sticky = Tk.W)

        self.DisplayedTypesFrame = Tk.Frame(self.MainWindow)
        self.DisplayedTypesFrame.grid(row = 3, column = 1, sticky = Tk.W)

        self.StreamsFrame = Tk.Frame(self.MainWindow)
        self.StreamsFrame.grid(row = 4, column = 1, columnspan = 2, sticky = Tk.W)

        self.StreamInfoFrame = Tk.Frame(self.MainWindow)
        self.StreamInfoFrame.grid(row = 0, column = 1, sticky = Tk.W)

        TauFrame = Tk.Frame(self.StreamInfoFrame)
        TauFrame.pack(anchor = Tk.W)

        # Initialize Event specific features
        self._EventDecryptFunctions = {}
        self._InitEventsVars()
        self._InitTrackerEventsVars()

        self.StreamsVariable = Tk.IntVar(master = self.MainWindow)
        self._MinimumGeometry = np.array([10, 10, 2])
        self.AddStream('None', self._MinimumGeometry)

        self.Display = Figure(figsize=(7,7), dpi=150)
        self.DisplayAx = self.Display.add_subplot(111)
        self.DisplayCMap = 'binary'
        self.DisplayImShow = self.DisplayAx.imshow(self.StreamsMaps['None'][:,:,0], vmin = 0, vmax = 1, origin = "lower", cmap = self.DisplayCMap)
        
        self.DisplayCanvas = FigureCanvasTkAgg(self.Display, self.DisplayFrame)
        self.DisplayCanvas.draw()
        self.DisplayCanvas.get_tk_widget().grid(row = 0, column = 0)
        
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

        self.RTDisplayCanvas = FigureCanvasTkAgg(self.RealTimeDisplay, self.RTFrame)
        self.RTDisplayCanvas.draw()
        self.RealTimeDisplay.tight_layout()
        self.RTDisplayCanvas.get_tk_widget().grid(row = 0, column = 0)

        self.FPSValue = 0.
        self.FPSLabel = Tk.Label(self.RTFrame, text = "FPS : {0}".format(0))
        self.FPSLabel.grid(row = 1, column = 0)
        
        self.ColorsSelection = ['Single color', 'Rainbow mode']
        self.ColorsVMax = [1., 2.]
        self.ColorsMode = 1
        self.ColorSelectionButton = Tk.Button(self.PolaritiesFrame, text = self.ColorsSelection[self.ColorsMode], command = self.ChangeColorMode)
        self.ColorSelectionButton.pack(anchor = Tk.W)

        L = Tk.Label(self.PolaritiesFrame, text = "Polarities selection :")
        L.pack(anchor = Tk.W)
        PolasModes = [("OFF", 0),
                     ("ON", 1)]
        self.DisplayedPolas = {Polarity: True for PolaName, Polarity in PolasModes}
        self.DisplayedPolaritiesVars = []

        self.PolaritiesButtons = []
        for nPola, PolaParams in enumerate(PolasModes):
            Text, Polarity = PolaParams
            self.DisplayedPolaritiesVars += [Tk.IntVar(master = self.MainWindow)]
            self.DisplayedPolaritiesVars[-1].set(int(self.DisplayedPolas[Polarity]))
            self.PolaritiesButtons += [Tk.Checkbutton(self.PolaritiesFrame, text = Text, variable = self.DisplayedPolaritiesVars[-1], command = lambda n=nPola, P=Polarity: self.SwitchDisplayedPolasOnOff(n,P))]
            self.PolaritiesButtons[-1].pack(anchor = Tk.W)

        L = Tk.Label(self.DisplayedTypesFrame, text = "Display selection :")
        L.pack(anchor = Tk.W)
        DisplayedTypes = [("Events", 0),
                     ("Trackers", 1)]
        self.UpdateEventTypeMethods = {0:self.UpdateEvents, 1:self.UpdateTrackerEvents}

        self.DisplayedTypes = {EventType: True for EventName, EventType in DisplayedTypes}
        self.DisplayedTypesVars = []

        self.DisplayedTypesButtons = []
        for nType, EventParams in enumerate(DisplayedTypes):
            Text, EventType = EventParams
            self.DisplayedTypesVars += [Tk.IntVar(master = self.MainWindow)]
            self.DisplayedTypesVars[-1].set(int(self.DisplayedTypes[EventType]))
            self.DisplayedTypesButtons += [Tk.Checkbutton(self.DisplayedTypesFrame, text = Text, variable = self.DisplayedTypesVars[-1], command = lambda n=nType, E=EventType: self.SwitchDisplayedTypeOnOff(n,E))]
            self.DisplayedTypesButtons[-1].pack(anchor = Tk.W)

        L = Tk.Label(self.StreamsFrame, text = "Streams selection :")
        L.pack(anchor = Tk.W)
        self.StreamsButtons = {'None': Tk.Radiobutton(self.StreamsFrame, text = self.StreamsInfos['None'], variable = self.StreamsVariable, value = 0, command = lambda:self.SwitchStream(0))}
        self.StreamsButtons['None'].pack(anchor = Tk.W)
        self.CurrentStream = self.Streams[self.StreamsVariable.get()]

        L = Tk.Label(self.StreamInfoFrame, text = "Current Timestamp :")
        L.pack(anchor = Tk.W)
        self.TSLabel = Tk.Label(self.StreamInfoFrame)
        self.TSLabel.pack(anchor = Tk.W)
        L = Tk.Label(self.StreamInfoFrame, text = "Pipeline efficiency :")
        L.pack(anchor = Tk.W)
        self.EfficiencyLabel = Tk.Label(self.StreamInfoFrame)
        self.EfficiencyLabel.pack(anchor = Tk.W)

        self.TauLabel = Tk.Label(TauFrame, text = "Tau : {0} ms".format(self.Tau))
        self.TauLabel.grid(column = 0, row = 0)
        TauMinusButton = Tk.Button(TauFrame, text = ' - ', command = lambda: self.ChangeTau(-1))
        TauMinusButton.grid(column = 1, row = 0)
        TauPlusButton = Tk.Button(TauFrame, text = ' + ', command = lambda: self.ChangeTau(+1))
        TauPlusButton.grid(column = 2, row = 0)
        self.dTauLabel = Tk.Label(TauFrame, text = "dTau : {0:.2f} ms".format(self.dTau))
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
        self.MainWindow.bind('<e>', lambda event: self.SetDisplayMode(0))
        self.MainWindow.bind('<s>', lambda event: self.SetDisplayMode(1))
        self.MainWindow.bind('<b>', lambda event: self.SetDisplayMode(2))
        self.MainWindow.bind('<Down>', lambda event: self.SwitchStream(+1))
        self.MainWindow.bind('<Up>', lambda event: self.SwitchStream(-1))
        self.MainWindow.bind('<c>', lambda event: self.ChangeColorMode())
        self.MainWindow.protocol('WM_DELETE_WINDOW', self._on_closing)

        ClearStreamsButton = Tk.Button(self.MainWindow, text='Clear', command = self.ClearStreams)
        ClearStreamsButton.grid(row = 5, column = 1, sticky = Tk.W)

        self.MainThread = MainThreadClass(self, self.MainSocket)
        self.MainThread.start()

        self.UpdateStreamsList()
        self.UpdateScene()
        self.UpdateTS()
        self.UpdateRTD()
        self.MainWindow.mainloop()

    def AddStream(self, Name, Geometry, StreamInfo = "No display", StreamType = 'Stream'):
        Geometry = np.maximum(Geometry, self._MinimumGeometry)
        self.Streams += [Name]
        self.StreamsTypes[Name] = StreamType
        self.StreamsTimes[Name] = 0
        self.LastTSAtRTDUpdate[Name] = 0
        self.StreamsInfos[Name] = StreamInfo
        self.AddEventsStreamVars(Name, Geometry)
        self.AddTrackersStreamVars(Name, Geometry)

    def RemoveStream(self, Name):
        del self.StreamsTimes[Name]
        del self.LastTSAtRTDUpdate[Name]
        del self.StreamsInfos[Name]
        del self.StreamsTypes[Name]
        self.RemoveEventsStreamVars(Name)
        self.RemoveTrackersStreamVars(Name)
        self.CurrentlyDisplayedStreams.remove(Name)
        self.Streams.remove(Name)
        self.CurrentStream = self.Streams[-1]

    def ChangeColorMode(self):
        self.ColorsMode = 1-self.ColorsMode
        self.ColorSelectionButton.configure(text = self.ColorsSelection[self.ColorsMode])
        self.DisplayImShow.set_clim(vmax = self.ColorsVMax[self.ColorsMode])

    def SwitchStream(self, var):
        if var != 0:
            self.StreamsVariable.set(min(max(self.StreamsVariable.get() + var, 0), len(self.Streams)-1))
        self.Reload()

    def SwitchDisplayedTypeOnOff(self, nType, EventType):
        if self.DisplayedTypesVars[nType].get():
            self.DisplayedTypes[EventType] = True
        else:
            self.DisplayedTypes[EventType] = False
        self.Reload()

    def SwitchDisplayedPolasOnOff(self, nPola, Polarity):
        if self.DisplayedPolaritiesVars[nPola].get():
            self.DisplayedPolas[Polarity] = True
        else:
            self.DisplayedPolas[Polarity] = False
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
        self.Tau = max(self.dTau, self.Tau + var * self.dTau)
        self.TauLabel['text'] = "Tau : {0} ms".format(self.Tau)

    def ChangedTau(self, var):
        self.dTau *= 2**var
        self.dTauLabel['text'] = "dTau : {0:.2f} ms".format(self.dTau)

    def printPolaritiesVariable(self):
        print(self.PolaritiesVariable.get())

    def Reload(self):
        self.CurrentStream = self.Streams[self.StreamsVariable.get()]
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
        self.CurrentRTValue = max(self.RTMinLim, 1000 * TimeEllapsed / self.UpdateRTDTimeConstant)

        LogValue = np.log10(self.CurrentRTValue)
        self.RTRectangle.set_color(self.RTDColors[max(min(int(10 * (LogValue + 4) / 5), 9), 0)])
        self.RTRectangle.set_height(self.CurrentRTValue)
        self.RealTimeDisplayAx.set_title("{0:.1f}".format(LogValue))

        self.RealTimeDisplay.canvas.draw()
        self.MainWindow.after(self.UpdateRTDTimeConstant, self.UpdateRTD)

    def UpdateScene(self):
        if self.ReloadMap:
            self.DisplayAx.clear()
        for EventType, IsActive in self.DisplayedTypes.items():
            self.UpdateEventTypeMethods[EventType](IsActive)

        self.ReloadMap = False

        t = time.time()
        self.EfficiencyLabel['text'] = "{0} ev/s".format(int(self.PreviousNEvents / (t - self.PreviousTime)))
        self.FPSValue = self.FPSValue * 0.8 + 0.2/(t - self.PreviousTime)
        self.FPSLabel['text'] = "FPS : {0}".format(int(self.FPSValue))
        self.PreviousTime = t
        self.PreviousNEvents = 0
        self.MainWindow.after(5, self.UpdateScene)

    def _DecryptFullEvent(self, Socket, eventList):
        t = eventList[1]
        self.StreamsTimes[Socket] = max(self.StreamsTimes[Socket], t)
        for Extension in eventList[2:]:
            self._EventDecryptFunctions[Extension[0]](Socket, t, Extension[1:])

    def _InitEventsVars(self):
        self.StreamsShapes = {}
        self.StreamsMaps = {}
        L = Tk.Label(self.EventSpecificFrame, text = "Currently displayed events :")
        L.pack(anchor = Tk.W)
        self.EventsLabel = Tk.Label(self.EventSpecificFrame)
        self.EventsLabel.pack(anchor = Tk.W)
        self._EventDecryptFunctions[1] = self._DecryptEvent

    def _DecryptEvent(self, Socket, t, eventList):
        self.StreamsMaps[Socket][eventList[0][0], eventList[0][1], eventList[1]] = t

    def AddEventsStreamVars(self, Name, Geometry):
        self.StreamsShapes[Name] = Geometry
        self.StreamsMaps[Name] = -10*np.ones(self.StreamsShapes[Name])
    def RemoveEventsStreamVars(self, Name):
        del self.StreamsMaps[Name]
        del self.StreamsShapes[Name]

    def UpdateEvents(self, ActiveType):
        if ActiveType:
            Map = np.zeros(self.StreamsShapes[self.CurrentStream][:2])
            StreamMap = self.StreamsMaps[self.CurrentStream]
            StreamTime = self.StreamsTimes[self.CurrentStream]
            for Polarity, IsActive in self.DisplayedPolas.items():
                if IsActive:
                    Map = (Map + ((StreamTime - StreamMap[:,:,Polarity]) < self.Tau/1000.)) #For color mode

            if self.ReloadMap:
                self.DisplayImShow = self.DisplayAx.imshow(np.transpose(Map), vmin = 0, vmax = self.ColorsVMax[self.ColorsMode], origin = "lower", cmap = self.DisplayCMap)
                self.DisplayAx.set_xlim(0, self.StreamsShapes[self.CurrentStream][0])
                self.DisplayAx.set_ylim(0, self.StreamsShapes[self.CurrentStream][1])
            else:
                self.DisplayImShow.set_data(np.transpose(Map))
            self.EventsLabel['text'] = "{0}".format(int(Map.sum()))
            self.Display.canvas.draw()
        elif self.ReloadMap:
            Map = np.zeros(self.StreamsShapes[self.CurrentStream][:2])
            self.EventsLabel['text'] = "0"
            Map[0,0] = 1
            self.DisplayImShow = self.DisplayAx.imshow(np.transpose(Map), vmin = 0, vmax = 1., origin = "lower", cmap = self.DisplayCMap)
            self.Display.canvas.draw()

    def _InitTrackerEventsVars(self):
        self.TrackersScalingSize = 20
        self.ActiveTrackers = {'None': {}}
        self.UpdatedTrackers = {}
        self.DisplayedTrackers = {}
        self.DisplayedTrackersIDs = {}
        self.CurrentStreamLastUpdatedAlphasTs = {}
        L = Tk.Label(self.EventSpecificFrame, text = "Active trackers :")
        L.pack(anchor = Tk.W)
        self.TrackersLabel = Tk.Label(self.EventSpecificFrame)
        self.TrackersLabel.pack(anchor = Tk.W)
        self._EventDecryptFunctions[2] = self._DecryptTrackerEvent

    def _DecryptTrackerEvent(self, Socket, t, eventList):
# Structure is [x, y, t, p, 1, ID, xTracker, yTracker, rTracker, sTracker, color, marker]
# Structure is [Location, ID, Angle, Scaling, Color, Marker]
        ID = eventList.pop(1)
        self.ActiveTrackers[Socket][ID] = eventList + [t,t]
        self.UpdatedTrackers[Socket].add(ID)

    def AddTrackersStreamVars(self, Name, Geometry):
        self.ActiveTrackers[Name] = {}
        self.UpdatedTrackers[Name] = set()
        self.CurrentStreamLastUpdatedAlphasTs[Name] = -np.inf
    def RemoveTrackersStreamVars(self, Name):
        del self.ActiveTrackers[Name]
        del self.UpdatedTrackers[Name]
        del self.CurrentStreamLastUpdatedAlphasTs[Name]
        #for TrackerID in list(self.DisplayedTrackers.keys()):
        #    self.RemoveTrackerParts(TrackerID)
    
    def AddTrackerParts(self, TrackerID, TrackerData, Alpha):
        self.DisplayedTrackers[TrackerID] = [None, None]
        self.DisplayedTrackers[TrackerID][0] = self.DisplayAx.plot(TrackerData[0][0], TrackerData[0][1], color = TrackerData[3], marker = TrackerData[4], alpha = Alpha)[0]
        Theta, S = TrackerData[1], TrackerData[2]
        dx, dy = self.TrackersScalingSize * S * np.cos(Theta), self.TrackersScalingSize * S * np.sin(Theta)
        self.DisplayedTrackers[TrackerID][1] = self.DisplayAx.plot([TrackerData[0][0], TrackerData[0][0]+dx], [TrackerData[0][1], TrackerData[0][1]+dy], color = TrackerData[3], alpha = Alpha)[0]
        self.DisplayedTrackersIDs[TrackerID] = self.DisplayAx.text(TrackerData[0][0] + 5, TrackerData[0][1], str(TrackerID), color = TrackerData[3], alpha = Alpha)

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

    def UpdateTrackerEvents(self, ActiveType):
        if ActiveType:
            T = self.StreamsTimes[self.CurrentStream]
            if self.ReloadMap:
                self.DisplayedTrackers = {}
                self.DisplayedTrackersIDs = {}
                for TrackerID, TrackerData in dict(self.ActiveTrackers[self.CurrentStream]).items():
                    Alpha = min(1, max(0, np.e**((TrackerData[-2] - T)/(self.Tau / 1000))))
                    self.AddTrackerParts(TrackerID, TrackerData, Alpha)
                self.UpdatedTrackers[self.CurrentStream].clear()
            else:
                UpdateAlphas = (T - self.CurrentStreamLastUpdatedAlphasTs[self.CurrentStream] > self.Tau / 5000)
                if UpdateAlphas:
                    self.CurrentStreamLastUpdatedAlphasTs[self.CurrentStream] = T
                    for TrackerID in list(self.DisplayedTrackers.keys()):
                        TrackerData = self.ActiveTrackers[self.CurrentStream][TrackerID]
                        Alpha = min(1, max(0, np.e**((TrackerData[-2] - T)/(self.Tau / 1000))))
                        if Alpha < 0.001:
                            self.RemoveTrackerParts(TrackerID)
                            continue
                        if TrackerID in self.UpdatedTrackers[self.CurrentStream]:
                            if TrackerID not in self.DisplayedTrackers.keys():
                                self.AddTrackerParts(TrackerID, TrackerData, Alpha)
                            else:
                                self.UpdateTrackerParts(TrackerID, TrackerData)
                            self.UpdateTrackerAlpha(TrackerID, Alpha)
                else:
                    for TrackerID in list(self.UpdatedTrackers[self.CurrentStream]):
                        TrackerData = self.ActiveTrackers[self.CurrentStream][TrackerID]
                        Alpha = min(1, max(0, np.e**((TrackerData[-2] - T)/(self.Tau / 1000))))
                        if TrackerID not in self.DisplayedTrackers.keys():
                            self.AddTrackerParts(TrackerID, TrackerData, Alpha)
                        else:
                            if Alpha < 0.001:
                                self.RemoveTrackerParts(TrackerID)
                            else:
                                self.UpdateTrackerParts(TrackerID, TrackerData)
                                self.UpdateTrackerAlpha(TrackerID, Alpha)
                self.UpdatedTrackers[self.CurrentStream].clear()
        self.TrackersLabel['text'] = "{0}".format(len(self.ActiveTrackers[self.CurrentStream]))

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
                if s in self.Parent.StreamsMaps.keys():
                    self.Parent.StreamsMaps[s] = -10*np.ones(self.Parent.StreamsShapes[s])
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
                self.Parent.StreamsMaps[Stream][self.Parent.StreamsMaps[Stream] > tNew] = -10
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
            data, addr = self.Connexion.recvfrom(2048)
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

