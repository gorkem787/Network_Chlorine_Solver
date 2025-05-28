# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 23:48:00 2024

@author: gorkem
"""

import numpy as np
import matplotlib.pyplot as plt
from epyt import epanet
import pandas as pd
from ADE_solver import solver
import timeit

start = timeit.default_timer()


# import warnings
# warnings.filterwarnings("ignore")


class node:
    def __init__(self, node_name, K, C0):
        self.node_name = node_name
        self.flow = 0
        self.n_calculation = 0
        self.K = K  # Bozunma katsayısı
        self.Co = C0
        self.consantration_w = np.array([C0])
        self.consantration_t = np.array([0])
        self.flow = np.array([0])
        self.node_type = "node"

    def set_consantration(self, consantration_t, consantration_w, flow):
        self.consantration_t = consantration_t  # Zaman dilimleri
        self.consantration_w = consantration_w  # Konsantrasyon değerleri
        self.flow = flow
        # self.n_calculation += 1

    def get_consantration(self):
        return self.consantration_t, self.consantration_w

    def add_consantration(self, t, c, Flow):
        if c == 0:
            print("Node Cons 0")
            print(self.node_name, c, round(t / 3600))
        index = np.flatnonzero(self.consantration_t == t)
        if self.node_type == "reservoir":
            pass
        elif self.node_type == "pump":
            pass
        else:
            if len(index) > 0:
                index = index[0]
                self.consantration_w[index] = (self.consantration_w[index] * self.flow[index] + c * Flow) / (
                            self.flow[index] + Flow)
                self.flow[index] += Flow
            else:
                self.consantration_t = np.append(self.consantration_t, np.array([t]))
                self.consantration_w = np.append(self.consantration_w, np.array([c]))
                self.flow = np.append(self.flow, np.array([Flow]))
        self.n_calculation += 1

    def add_reservoir(self, initial_volume):
        self.node_type = "reservoir"
        self.initial_volume = initial_volume

    def add_pump(self):
        self.node_type = "pump"

    def plot(self):  # Konsantrasyonun zamanla değişimini grafikte gösterir
        fig, ax = plt.subplots()
        # SMALL_SIZE = 5
        # MEDIUM_SIZE = 5
        # BIGGER_SIZE = 7

        # plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        # plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        # plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        # plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        # plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        # plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        # plt.rc('figure', titlesize=BIGGER_SIZE)

        ax.plot(self.consantration_t, self.consantration_w)
        ax.set_xlabel('Time (seconds)')
        ax.set_ylabel('Concentration(g/m^3)')
        ax.set_title(f'Concentration Over Time at node{self.node_name}')
        plt.show()


# Pipe sınıfı boru özelliklerini tutar ve konsantrasyon hesaplarını yapar

class pipe():
    def __init__(self, number, flow, length, diameter, nodes, t_step, part_lenght):  # id, m3/s, m, m, nodes
        self.s_node = nodes[0]  # Başlangıç düğümü
        self.d_node = nodes[1]
        self.flow = np.array(flow)
        self.length = length  # Borunun uzunluğu
        self.U_list = 4 * self.flow / (np.pi * diameter ** 2)  # Akış hızı
        if min(abs(self.U_list)) <= 0.0003:
            print("Pipe ıd slow " + str(number))
            print(self.U_list)

        self.number = number
        self.diameter = diameter
        self.t_step = t_step
        self.make_parts(part_lenght)

    def t_step(self, t_step):
        self.t_step = t_step

    def t_index(self, t):
        t = round(t / self.t_step)
        return t

    def make_parts(self, part_lenght):
        self.part_lenght = part_lenght

        self.sub_node_count = round(self.length / self.part_lenght)
        self.part_lenght = self.length / self.sub_node_count
        # print(self.part_lenght)
        self.subnode_list = []
        self.subnode_list.append(self.s_node)
        for n in range(self.sub_node_count):
            sub_node = node(f"sub_node{n}", K=self.s_node.K, C0=self.s_node.Co)
            # sub_node.add_consantration(0, self.s_node.Co, 1e5)
            self.subnode_list.append(sub_node)
        self.subnode_list.append(self.d_node)
        return self.part_lenght

    def check_direction(self, t):
        U = self.U_list[self.t_index(t)]
        # print(U)
        if U < 0:
            temp = self.s_node
            self.s_node = self.d_node
            self.d_node = temp
            self.U_list = -1 * self.U_list
            self.flow = -1 * self.flow
            self.subnode_list.reverse()
        return self.s_node, self.d_node

    def calculate_sub_nodes(self, t, flow, U, x, x_s, Co, diameter, K):
        for sub_node_index in range(len(self.subnode_list) - 1):
            index = self.t_index(t)
            # print(index)
            source_t, source_ws = self.subnode_list[sub_node_index].get_consantration()
            source_t = source_t[[index - 1, index]]
            source_ws = source_ws[[index - 1, index]]
            _, Co = self.subnode_list[sub_node_index].get_consantration()
            Co = Co[index - 1]
            if U > 3e-200:
                c_t = solver(t=t, U=U, x=x, x_source=x_s, Co=Co, source_t=source_t,
                             source_ws=source_ws, diameter=diameter, K=K)
                self.subnode_list[sub_node_index + 1].add_consantration(t, c_t, flow)
            elif U < 3e-200:
                source_t, source_ws = self.subnode_list[sub_node_index + 1].get_consantration()
                c_t = (1 - K) * source_ws[[index - 1]]
                # print(c_t)
                self.subnode_list[sub_node_index + 1].add_consantration(t, c_t, flow)

    def calculate_time(self, t):
        flow = self.flow[self.t_index(t)]
        U = self.U_list[self.t_index(t)]
        x = self.part_lenght
        x_source = 0
        Co = self.s_node.Co
        K = self.s_node.K
        diameter = self.diameter
        self.calculate_sub_nodes(t=t, flow=flow, U=U, x=x, x_s=x_source, Co=Co, diameter=diameter, K=K)


class system(node, pipe):
    def __init__(self, max_time, time_step, pipes=[]):
        self.pipes = pipes
        self.max_time = max_time
        self.time_step = time_step

    def pipe_connections(self, t):
        connetion_data = []
        for p in self.pipes:
            s_node, d_node = p.check_direction(t)
            connetion_data.append([p, s_node, d_node])
        connetion_data = np.array(connetion_data)
        return connetion_data

    def system_calculate_time(self):
        time_step = self.time_step
        max_time = self.max_time
        t = time_step
        while t <= max_time:
            print(t)
            connection_data = self.pipe_connections(t)
            pipe_calculated = np.zeros(len(connection_data[:, 0]))
            for node in connection_data[:, 1:2].flat:
                node.n_calculation = 0
            while pipe_calculated.all() != 1:
                for index in range(len(connection_data[:, 0])):
                    n_pipe = connection_data[index, 0]
                    n_source = connection_data[index, 1]
                    n_destination = connection_data[index, 2]
                    con_pipe = np.sum(np.equal(connection_data[:, 2], n_source))
                    non_calculated_pipe = con_pipe - n_source.n_calculation
                    if non_calculated_pipe == 0:
                        if pipe_calculated[index] == 0:
                            n_pipe.calculate_time(t)
                            pipe_calculated[index] = 1
            t += time_step

    def constructor(self, file_name, K=1 / 86400, part_lenght=1):
        # Initialize EPANET with the provided input file    
        self.K = K
        self.file_name = file_name
        net = epanet(file_name)
        units = net.getUnits().to_dict()
        net.setTimeReportingStep(self.time_step)
        # code only valid for SI units
        if units["Units_US_Customary"] == 1:
            GPM_m3s = 6.309E-5
            inc_m = 0.0254
            feet_m = 0.3048
        else:
            GPM_m3s = 0.001
            inc_m = 0.001
            feet_m = 1
        # Get network information   
        pipe_links = net.getNodesConnectingLinksIndex()
        node_demand = (net.getNodeBaseDemands())[1]
        pump_node = [i for i, demand in enumerate(node_demand) if demand < 0][0] + 1
        reservoir_node = net.getNodeTankIndex()[0]

        print(reservoir_node)
        print(pump_node)
        # works only for 1 pumps
        nodes_index = net.getNodeIndex()
        pipe_diameter = net.getLinkDiameter()
        pipe_len = net.getLinkLength()
        initial_quality = np.array(net.getNodeInitialQuality())

        net.openHydraulicAnalysis()
        net.openQualityAnalysis()
        net.initializeHydraulicAnalysis(0)
        net.initializeQualityAnalysis(0)
        tstep = 1
        Time, Flows, Quality_pump, reservoir = [], [], [], []
        while tstep > 0:
            time = net.runHydraulicAnalysis()
            quality = net.runQualityAnalysis()
            Flows.append(net.getLinkFlows())
            Quality_pump.append(net.getNodeActualQuality(pump_node))
            reservoir.append(net.getNodeActualQuality(reservoir_node))
            Time.append(time)
            tstep = net.nextHydraulicAnalysisStep()
            qtstep = net.nextQualityAnalysisStep()

        # Plot the network (optional)
        net.plot(nodesID=True)
        net.closeQualityAnalysis()
        net.closeHydraulicAnalysis()

        # print(len(initial_quality))
        Time = np.array(Time)
        Flows = np.array(Flows)
        Quality_pump = np.array(Quality_pump)
        reservoir = np.array(reservoir)

        # Create nodes for the network
        self.node_list = []
        for i in range(len(nodes_index)):
            i += 1
            if i == pump_node:
                n = node(i, K, initial_quality[i - 1])
                # print(Quality_pump)
                n.set_consantration(np.array(Time), Quality_pump, flow=np.zeros(len(Time)))
                n.add_pump()
            elif i == reservoir_node:
                n = node(i, K, initial_quality[i - 1])
                n.set_consantration(np.array(Time), reservoir, flow=np.zeros(len(Time)))
                n.add_reservoir(1)
            else:
                n = node(i, K, initial_quality[i - 1])
            self.node_list.append(n)

            # Create pipes for the network
        self.pipes = []
        for i in range(len(pipe_links)):
            flow = Flows[:, i] * GPM_m3s
            lenght = pipe_len[i] * feet_m
            diameter = pipe_diameter[i] * inc_m
            source_node_id = pipe_links[i][0]
            destination_node_id = pipe_links[i][1]
            source_node = self.node_list[int(source_node_id) - 1]
            destination_node = self.node_list[int(destination_node_id) - 1]
            # print((i+1, flow,lenght, diameter, [source_node.node_name,destination_node.node_name]))
            self.pipes.append(pipe(i + 1, flow, lenght, diameter,
                                   [source_node, destination_node], t_step=self.time_step, part_lenght=part_lenght))
        return self.node_list, self.pipes

    def plot_with_epanet(self, save=None):

        net = epanet(self.file_name)
        net.openHydraulicAnalysis()
        net.openQualityAnalysis()
        net.setTimeReportingStep(self.time_step)
        net.initializeHydraulicAnalysis(0)
        net.initializeQualityAnalysis(0)
        tstep = 1
        Time, Flows, quality = [], [], []
        while tstep > 0:
            time = net.runHydraulicAnalysis()
            q1 = net.runQualityAnalysis()
            Flows.append(net.getLinkFlows())
            quality.append(net.getNodeActualQuality())
            Time.append(time)
            tstep = net.nextHydraulicAnalysisStep()
            qtstep = net.nextQualityAnalysisStep()
            self.save_excel = []
        index = 0
        Time = np.array(Time)
        for n in self.node_list:
            c2 = np.array(quality)[:, index]
            t, c1 = n.get_consantration()
            fig, ax = plt.subplots()
            ax.plot(t / 3600, c1, label="mazaheri")
            ax.plot(Time / 3600, c2, label="epanet")
            ax.set_xlabel('Time (hour)')
            ax.set_ylabel('Concentration(g/m^3)')
            ax.set_title(f'Concentration Over Time at {n.node_type, net.getNodeNameID(n.node_name)} K={self.K}')
            ax.legend()
            if save == True:
                self.save_excel.append([net.getNodeNameID(n.node_name), t / self.time_step, c1, c2])
            index += 1
            fig.show()

a = system(max_time=55 * 3600, time_step=3600 / 4)
nodes, pipes = a.constructor("epanet_net2/net2.inp", K=1.5 / 86400)
# nodes, pipes = a.constructor("h_deneme/h_deneme.inp", K = 0/86400)
# nodes, pipes = a.constructor("y_deneme/y_deneme.inp", K = 0/86400)
a.system_calculate_time()
print("The difference of time is :",
      timeit.default_timer() - start)
a.plot_with_epanet(save=False)
# df = pd.DataFrame(data = a.save_excel, columns = ["node_id","Time", "Mazaheri", "Epanet"]).explode(["Time", "Mazaheri", "Epanet"])
# df.index += 1
# df.to_excel("mazaheri_degerler_k=1.5.xlsx")
# for pipe in p:
#     print([pipe.s_node.node_name,pipe.d_node.node_name])
# for node in n:
# node.plot()
