import numpy as np
import socket
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import _pickle as cPickle
import tkinter as Tk
from tkinter import filedialog as tkFileDialog
import time
from matplotlib.patches import Rectangle
from colour import Color
import random

import threading
import datetime
from string import ascii_lowercase as alphabet

import EventsHandlers

matplotlib.use("TkAgg")

class Display:
    _DefaultTau = 0.030
    def __init__(self, mainPort = 54242, questionPort = 54243, responsePort = 54244, listen = ""):
        self.SharedData = {'t':0, 'Tau':self._DefaultTau, 'Stream':None, 'ReloadCommand': self.Reload}
        self.dTau = 0.005
        self.PreviousTime = 0.
        self.PreviousSceneTime = 0.
        self.PreviousNEvents = 0
	
        self.MainSocket = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
        self.MainSocket.bind((listen, mainPort))

        self.QuestionSocket = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
        self.QuestionSocket.bind((listen, questionPort))

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
        for nType, Handler in enumerate(EventsHandlers._HANDLERS):
            if Handler.DisplayType != 'Info':
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
            elif Handler.DisplayType == 'Info':
                HandlerInfoFrame = self.HandlersMapInfoFrame
                HandlerFrame = None
                Button = None
            else:
                raise Exception("Ill defined handler")
            if not HandlerFrame is None:
                Button.grid(row = HandlersRows[Handler.DisplayType], column = 0, sticky = Tk.W)
                HandlerOptionsFrame = Tk.Frame(HandlerFrame)
                HandlerOptionsFrame.grid(row = HandlersRows[Handler.DisplayType], column = 1, sticky = Tk.W)
                HandlersRows[Handler.DisplayType] += 1
            else:
                HandlerOptionsFrame = None

            HandlerInfoFrame = Tk.Frame(HandlerInfoFrame)
            HandlerInfoFrame.pack(anchor = Tk.W)
            
            self.Handlers[Handler.Key] = Handler(self.MainWindow, self.DisplayAx, self.Display.canvas, HandlerInfoFrame, HandlerOptionsFrame, self.SharedData)
        for UnusedKey in range(10):
            if UnusedKey not in self.Handlers:
                self.Handlers[UnusedKey] = EventsHandlers.NullHandler()

        self.StreamsVariable = Tk.IntVar(master = self.MainWindow)
        self._MinimumGeometry = np.array([10, 10])
        self.AddStream('None', self._MinimumGeometry)

        self.RTMinLim = 10**-4
        self.RTMaxLim = 10**1
        self.CurrentRTValue = 1
        self.RealTimeDisplay = Figure(figsize=(1,2), dpi=100)
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

        self.TSLabel = Tk.Label(self.StreamInfoFrame, text = "Current timestamp : ")
        self.TSLabel.pack(anchor = Tk.W)
        self.EfficiencyLabel = Tk.Label(self.StreamInfoFrame, text = "Pipeline efficiency :")
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

        self._DelayedDataGCTargets = []

        self.QuestionThread = QuestionThreadClass(self, self.QuestionSocket, socket.socket(socket.AF_INET,socket.SOCK_DGRAM), responsePort)
        self.QuestionThread.start()

        self.MainThread = MainThreadClass(self, self.MainSocket)
        self.MainThread.start()

        self.UpdateRTDTimeConstant = 1000 # in ms
        self.UpdateSceneTimeConstant = 5
        self.UpdateTSTimeConstant = 50
        self.UpdateStreamsListTimeConstant = 100
        self.FullUpdateCheckTimeConstant = 3000

        self.LastLoopUpdates = {'StreamsList':[0,0], 'Scene':[0,0], 'TS':[0,0], 'RTD':[0,0]}
        self.FullUpdateCheck()
        self.MainWindow.mainloop()

    def FullUpdateCheck(self):
        for Loop, (LastUpdateID, LastRegisteredUpdateID) in self.LastLoopUpdates.items():
            if LastUpdateID != LastRegisteredUpdateID: # An update occured, meaning that loop is still running
                self.LastLoopUpdates[Loop][1] = self.LastLoopUpdates[Loop][0] # We aknowledge that update
            else:
                print("Re-starting loop {0}".format(Loop))
                try:
                    self.__class__.__dict__['LoopUpdate'+Loop](self) # Else we restart that loop, and do not aknoledge any update. Next occurence, LastUpdateID should have changed
                except Exception as e:
                    print("Restart failure", e)
        self.MainWindow.after(self.FullUpdateCheckTimeConstant, self.FullUpdateCheck)

    def LoopUpdateStreamsList(self, Cycle = True):
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
        if Cycle:
            self.LastLoopUpdates['StreamsList'][0] = random.random()
            self.MainWindow.after(self.UpdateStreamsListTimeConstant, self.LoopUpdateStreamsList)

    def LoopUpdateRTD(self):
        TimeEllapsed = self.StreamsTimes[self.CurrentStream] - self.LastTSAtRTDUpdate[self.CurrentStream]
        self.LastTSAtRTDUpdate[self.CurrentStream] = self.StreamsTimes[self.CurrentStream]
        self.CurrentRTValue = max(self.RTMinLim, TimeEllapsed / self.UpdateRTDTimeConstant * 1e3)

        LogValue = np.log10(self.CurrentRTValue)
        self.RTRectangle.set_color(self.RTDColors[max(min(int(10 * (LogValue + 4) / 5), 9), 0)])
        self.RTRectangle.set_height(self.CurrentRTValue)
        self.RealTimeDisplayAx.set_title("{0:.1f}".format(LogValue))

        self.RealTimeDisplay.canvas.draw()
        self.LastLoopUpdates['RTD'][0] = random.random()
        t = time.time()
        self.EfficiencyLabel['text'] = "Pipeline efficiency : {0} ev/s".format(int(self.PreviousNEvents / (t - self.PreviousTime)))
        self.PreviousTime = t
        self.PreviousNEvents = 0
        self.MainWindow.after(self.UpdateRTDTimeConstant, self.LoopUpdateRTD)

    def LoopUpdateScene(self):
        if self.ReloadMap:
            self.DisplayAx.clear()
        self.SharedData['t'] = self.StreamsTimes[self.CurrentStream]
        for Handler in self.Handlers.values():
            Handler.Update(self.ReloadMap)

        self.Display.canvas.draw()
        self.ReloadMap = False

        t = time.time()
        self.FPSValue = self.FPSValue * 0.9 + 0.1/(t - self.PreviousSceneTime)
        self.FPSLabel['text'] = "FPS : {0}".format(int(self.FPSValue))
        self.PreviousSceneTime = t
        self.LastLoopUpdates['Scene'][0] = random.random()
        self.MainWindow.after(self.UpdateSceneTimeConstant, self.LoopUpdateScene)

    def LoopUpdateTS(self):
        self.TSLabel['text'] = "Current timestamp : {0} ms".format(int(1000*self.StreamsTimes[self.CurrentStream]))
        self.LastLoopUpdates['TS'][0] = random.random()
        self.MainWindow.after(self.UpdateTSTimeConstant, self.LoopUpdateTS)

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
        if Name not in self.Streams:
            return
        self._DelayedDataGCTargets += [Name]
        for Stream in reversed(self.Streams):
            if Stream not in self._DelayedDataGCTargets:
                self.SharedData['CurrentStream'] = Stream
                break
        self.MainWindow.after(300, self.DelayedDataGC)

    def DelayedDataGC(self):
        Name = self._DelayedDataGCTargets.pop(0)
        del self.StreamsTimes[Name]
        del self.LastTSAtRTDUpdate[Name]
        del self.StreamsInfos[Name]
        del self.StreamsTypes[Name]
        del self.StreamsTaus[Name]
        for Handler in self.Handlers.values():
            Handler.DelStreamVars(Name)
        self.CurrentlyDisplayedStreams.remove(Name)
        self.Streams.remove(Name)
        self.ReloadNames = True

    def OnClick(self, event):
        if event.inaxes is None:
            return
        x, y = int(event.xdata), int(event.ydata)
        try:
            V = self.Handlers[self.DisplayedMapTypeVar.get()].DisplayImShow.get_array()[y, x] # inverted due to transpose
            self.PointLabel['text'] = "x = {0}, y = {1}, V = {2:.3f}".format(x, y, V)
        except Exception as e:
            print(e)
            self.PointLabel['text'] = "x = {0}, y = {1}, V = ?".format(x, y)

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
        self.UpdateStreamsList(Cycle = False)

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

        self.MainWindow.quit()
        self.MainWindow.destroy()

    def _DecryptFullEvent(self, Socket, eventDict):
        t = eventDict[0]
        self.StreamsTimes[Socket] = max(self.StreamsTimes[Socket], t)
        if 1 in eventDict:
            Location = eventDict[1][0]
        else:
            Location = None
        for ExtensionKey, Data in eventDict.items():
            if ExtensionKey == 0:
                continue
            self.Handlers[ExtensionKey].Decrypt(Socket, t, Data, Location)

class QuestionThreadClass(threading.Thread):
    def __init__(self, parent, QuestionConnexion, ResponseConnexion, ResponsePort):
        self.QuestionConnexion = QuestionConnexion
        self.ResponseConnexion = ResponseConnexion
        self.ResponsePort = ResponsePort
        self.Parent = parent
        threading.Thread.__init__(self)
        self.Run = True

    def run(self):
        while self.Run:
            data_raw, (addr, socket) = self.QuestionConnexion.recvfrom(1024)
            
            data = cPickle.loads(data_raw)
            print("Received {0} command from {1}, socket {2}".format(data['command'], addr, socket))
            print(data)
            if data['command'] == 'kill':
                self.Run = False
                continue
            Response = None
            Params = None
            if data['command'] == 'isup':
                Response = True
            if data['command'] == 'asksocket':
                idgiven = str(socket)
                if idgiven not in self.Parent.Streams:
                    Response = True
                    Params = idgiven
                    self.Parent.AddStream(idgiven, data['shape'], "Unknown ({0})".format(idgiven))
                else:
                    Response = False

            if data['command'] == 'askspecificsocket':
                socketasked = data['socket']
                Response = True
                Params = socketasked
                if socketasked not in self.Parent.Streams:
                    self.Parent.AddStream(socketasked, data['shape'], "Unknown ({0})".format(socketasked))

            if data['command'] == "cleansocket":
                s = data['socket']
                if s in self.Parent.StreamsTimes.keys():
                    self.Parent.StreamsTimes[s] = -0
                    self.Parent.LastTSAtRTDUpdate[s] = -0
                    Response = True
                else:
                    Response = False

            if data['command'] == 'socketdata':
                s = data['socket']
                if s in self.Parent.Streams:
                    self.Parent.StreamsInfos[s] = str(addr) + '@' + str(s) + ' : ' + data['infosline1']+"\n"+data['infosline2']
                    self.Parent.ReloadNames = True
                    Response = True
                else:
                    Response = False

            if data['command'] == 'destroysocket':
                s = data['socket']
                try:
                    self.Parent.RemoveStream(s)
                    Response = True
                except:
                    Response = False
            if Response is None:
                print("Question not answered")
            self.Respond(data['id'], (Response, Params), addr)

        self.QuestionConnexion.close()
        print("Questions connexion closed")

    def Respond(self, QuestionID, Response, addr):
        print("Answer to question {0} for {1}@{2} : {3}".format(QuestionID, addr, self.ResponsePort, Response))
        self.ResponseConnexion.sendto(cPickle.dumps({'id':QuestionID, 'answer':Response}), (addr, self.ResponsePort))

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

D = Display()

