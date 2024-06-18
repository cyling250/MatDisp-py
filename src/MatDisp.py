# -*- coding: utf-8 -*-
"""
    矩阵位移法计算核心模块
    调用顺序:
    1 初始化 Structure()
    2 调用 Node() 生成节点
    3 调用 Structure.link 生成杆件, 并记得保存后续可能会用到的杆件的引用

    4 调用 Element.add_strain() 为杆件添加约束
    5 调用 Element.set_load() 或 Node.set_load() 添加杆件荷载或节点荷载
    步骤4和步骤5没有强制的先后顺序

    此时，建模已经完成
    6 调用 Structure.get_freedom() 计算自由度
    7 调用 Structure.get_entire_k() 生成整体刚度矩阵
    8 调用 Structure.load_process() 生成荷载列阵
    9 调用 Structure.freedom_process() 缩减计算行列式
    此时，计算前准备已经完成
    10 使用 Structure.solve() 求解行列式
    11 使用 Structure.get_internal_force() 得到杆端内力

    @Author  : RichardoGu
    @Time    : 2024/6/18 16:31
"""
import numpy as np


class Node:
    """节点类"""
    index = 0  # 索引
    count = 0  # 个数

    def __init__(self,
                 x: float = 0,
                 y: float = 0):
        self.x = x
        self.y = y
        self.freedom = np.zeros(3, dtype=int)  # 节点自由度编号，一般有三个元素，半铰接点有特殊处理函数
        self.load = np.zeros(3)  # 荷载列阵,节点上的荷载只考虑前三个，节点平动和节点转动

        self.index = Node.index  # 节点编号
        Node.count += 1
        Node.index += 1

    def __eq__(self, other):
        """定义两个节点的相等"""
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        """实现Set"""
        return hash(self.index)

    def __str__(self):
        return f"[{self.x},{self.y}]"

    def __del__(self):
        Node.count -= 1

    def set_load(self, load: np.ndarray):
        """荷载列阵"""
        self.load += load


class Element:
    """单元类"""
    index = 0
    count = 0

    def __init__(self,
                 node1: Node = None,
                 node2: Node = None,
                 link_way1: bool = False,
                 link_way2: bool = False,
                 E: float = 1,
                 I: float = 1,
                 A: float = 1):
        """
        单元初始化函数
        :param node1: 节点1
        :param node2: 节点2
        :param link_way1: 节点1的连接方式
        :param link_way2: 节点2的连接方式
        :param E: 弹性模量
        :param I: 转动惯量
        :param A: 截面积
        """
        self.nodes: [Node, Node] = [node1, node2]  # 单元所包含的节点
        self.link_way: [bool, bool] = [link_way1, link_way2]  # 单元节点的连接方式 True为铰接,False为刚接
        self.freedom = np.zeros(6, dtype=int)  # 定位向量及自由度
        self.E: float = E  # 单元弹性模量
        self.I: float = I  # 单元截面惯性矩
        self.A = A  # 单元截面面积

        cosa = (self.nodes[1].x - self.nodes[0].x) / self.__len__()
        sina = (self.nodes[1].y - self.nodes[0].y) / self.__len__()
        self.T = np.array([
            [cosa, sina, 0, 0, 0, 0],
            [-sina, cosa, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, cosa, sina, 0],
            [0, 0, 0, -sina, cosa, 0],
            [0, 0, 0, 0, 0, 1]], dtype=float)  # 坐标变换矩阵

        self.k_e = self.T.T @ self.get_local_k_e() @ self.T  # 整体坐标系下单元刚度矩阵

        self.restrain = np.zeros(6, dtype=bool)  # 约束条件, False表示未约束，True表示约束
        self.load = np.array([0, 0, 0, 0, 0, 0], dtype=float)  # 局部坐标荷载列阵
        self.delta = np.zeros(6, dtype=float)  # 杆端位移
        self.force = np.zeros(6, dtype=float)  # 杆端内力

        self.index = Element.index  # 单元编号
        Element.count += 1  # 单元个数加1
        Element.index += 1  # 索引加1

    def __del__(self):
        Element.count -= 1

    def __eq__(self, other):
        """定义两个杆件相等"""
        if self.nodes[0] == other.nodes[0] and self.nodes[1] == other.nodes[1]:
            return True
        else:
            return False

    def __len__(self):
        """杆件单元的长度"""
        return ((self.nodes[0].x - self.nodes[1].x) ** 2 + (self.nodes[0].y - self.nodes[1].y) ** 2) ** 0.5

    def __hash__(self):
        return hash(self.index)

    def __str__(self):
        return f"{self.nodes[0].__str__()}->{self.nodes[1].__str__()}"

    def set_load(self, load):
        """
        转化等效节点力
        Load包含了各种荷载形式，荷载形式可以自行定义，只需要最终作用到self.Load变量上即可
        """
        length = self.__len__()
        if load[0] == 0:
            # 0表示杆上集中荷载，杆上集中荷载的形式[0,Fp,a],Fp是力的大小，a是集中力到杆件正方向左端的距离
            b = length - load[2]
            a = load[2]
            fp = load[1]
            self.load += np.array(
                [0, -fp * b ** 2 * (1 + 2 * a / length) / length ** 2, -fp * a * b ** 2 / length ** 2,
                 0, -fp * a ** 2 * (1 + 2 * b / length) / length ** 2, fp * a * b ** 2 / length ** 2], dtype=float)
        elif load[0] == 1:
            # 1表示杆端弯矩，Load[1]存放杆左端弯矩,Load[2]存放杆右端弯矩
            self.load[2] -= load[1]
            self.load[5] -= load[2]
        elif load[0] == 2:
            # 2表示均布荷载，Load[1]存放均布荷载的大小
            q = load[1]
            self.load += np.array([0., -q * length / 2., -q * length ** 2 / 12.,
                                   0., -q * length / 2., q * length ** 2 / 12.])
            # print("最初加上的荷载列阵为:", self.Load)
        elif load[0] == 3:
            # 3表示具有等斜度的均布荷载
            q = load[1]
            cosa = self.T[0][0]
            self.load += np.matmul(
                self.T,
                np.array([0., -q * length / 2. * cosa, -q * length ** 2 / 12. * cosa ** 2,
                          0., -q * length / 2. * cosa, q * length ** 2 / 12. * cosa ** 2]))

        return

    def get_local_k_e(self):
        E, A, L, I = self.E, self.A, self.__len__(), self.I
        if self.link_way[0] and self.link_way[1]:
            return np.array([
                [E * A / L, 0, 0, -E * A / L, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [-E * A / L, 0, 0, E * A / L, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0]
            ])
        else:
            return np.array([
                [self.E * A / L, 0, 0, -E * A / L, 0, 0],
                [0, 12 * E * I / pow(L, 3), 6 * E * I / pow(L, 2), 0, -12 * E * I / pow(L, 3), 6 * E * I / pow(L, 2)],
                [0, 6 * E * I / pow(L, 2), 4 * E * I / L, 0, -6 * E * I / pow(L, 2), 2 * E * I / L],
                [-E * A / L, 0, 0, E * A / L, 0, 0],
                [0, -12 * E * I / pow(L, 3), -6 * E * I / pow(L, 2), 0, 12 * E * I / pow(L, 3), -6 * E * I / pow(L, 2)],
                [0, 6 * E * I / pow(L, 2), 2 * E * I / L, 0, -6 * E * I / pow(L, 2), 4 * E * I / L]
            ])

    def add_restrain(self, restrain: [bool]):
        """为单元施加约束"""
        self.restrain = np.array(restrain)


class Structure:
    """结构类"""

    def __init__(self):
        self.elements: set[Element] = set()  # 单元集合
        self.size_of_K = 0  # 整体刚度矩阵维度
        self.in_result = None  # 位移定位矩阵，协助内力计算
        self.K = None  # 整体刚度矩阵
        self.load = None  # 荷载列阵
        self.result = None  # 结果矩阵

    def link(self,
             n1: Node,
             n2: Node,
             link_way1: bool = False,
             link_way2: bool = False,
             E: float = 1,
             I: float = 1,
             A: float = 1):
        """连接函数"""
        element = Element(n1, n2, link_way1, link_way2, E, I, A)
        self.elements.add(element)
        return element

    def get_link_dict(self) -> dict[Node, set[Element]]:
        """获得当前的连接状态"""
        link_dict = {}
        for i in self.elements:
            if i.nodes[0] not in link_dict:
                link_dict[i.nodes[0]] = set()
            if i.nodes[1] not in link_dict:
                link_dict[i.nodes[1]] = set()
            link_dict[i.nodes[0]].add(i)
            link_dict[i.nodes[1]].add(i)
        return link_dict

    def get_freedom(self):
        # 为每个杆件生成定位向量，并编号
        freedom_num = 0  # 初始化自由度编号
        link_dict: dict[Node, set[Element]] = self.get_link_dict()
        for node in link_dict.keys():
            link_temp = 0 - is_splice(node, link_dict[node])
            for element in link_dict[node]:
                i = element.nodes[1] == node  # 如果是第一个节点 则i=0,如果是第二个节点 则i=1
                if element.link_way[i]:
                    # 如果该节点与该单元的链接方式为铰接, 则其弯矩为 freedom_num + 2 + link_temp
                    link_temp += 1
                    element.freedom[i * 3 + 0] = freedom_num
                    element.freedom[i * 3 + 1] = freedom_num + 1
                    element.freedom[i * 3 + 2] = freedom_num + 2 + link_temp
                else:  # 如果连接方式为刚接, 则其弯矩为 freedom_num + 2
                    element.freedom[i * 3 + 0] = freedom_num
                    element.freedom[i * 3 + 1] = freedom_num + 1
                    element.freedom[i * 3 + 2] = freedom_num + 2

            if link_temp > 0:
                # 填充因为铰接而增加的转动自由度
                node.freedom = np.concatenate((node.freedom, np.zeros(link_temp, dtype=int)))
            for k in range(freedom_num, freedom_num + 3 + link_temp):
                node.freedom[k - freedom_num] = k
            freedom_num = freedom_num + 3 + link_temp
        return freedom_num

    def get_entire_k(self, force_calc: bool = False):
        """
        生成并获得整体刚度矩阵
        :param force_calc: 是否强制计算
        :return: 整体刚度矩阵
        """
        if self.K and not force_calc:
            return self.K

        self.K = np.zeros((self.size_of_K, self.size_of_K), dtype=float)
        for element in self.elements:
            for i in range(6):
                for j in range(6):
                    self.K[element.freedom[i]][element.freedom[j]] += element.k_e[i][j]  # 对号入座
        return self.K

    def load_process(self):
        """按照施加的荷载 生成整体荷载列阵"""
        self.load = np.zeros(self.size_of_K, dtype=float)
        # 先处理节点荷载
        for node in self.get_link_dict().keys():
            for j in range(3):
                self.load[node.freedom[0:3][j]] += node.load[j]  # 结点不需要反号
        # 再处理等效荷载
        for element in self.elements:
            element.load = element.T.T @ element.load  # 坐标转换，局部坐标转整体坐标
            for j in range(6):
                self.load[element.freedom[j]] += element.load[j]

    def freedom_process(self):
        """
        调用前置条件: 整体刚度矩阵已经形成 约束已经添加
        对整体刚度矩阵和荷载列阵进行划行划列
        :return:
        """
        # 将节点的约束添加到单元上
        delete_temp = []
        self.in_result = np.ones(self.size_of_K, dtype=bool)  # 初始化位移定位矩阵
        for element in self.elements:
            for j in range(6):
                if element.restrain[j]:  # 该自由度被约束住了
                    self.in_result[element.freedom[j]] = False
                    delete_temp.append(element.freedom[j])
            if element.link_way[0] and element.link_way[1]:
                # 说明这是一个轴力杆,对该轴力杆的弯矩进行自行约束
                self.in_result[element.freedom[2]] = False
                self.in_result[element.freedom[5]] = False
                delete_temp.append(element.freedom[2])
                delete_temp.append(element.freedom[5])
        if len(self.K) > 0:
            self.K = np.delete(self.K, [] if delete_temp is None else delete_temp, axis=0)  # 约束划行
            self.K = np.delete(self.K, [] if delete_temp is None else delete_temp, axis=1)  # 约束划列
        if len(self.load) > 0:
            self.load = np.delete(self.load, [] if delete_temp is None else delete_temp, axis=0)  # 约束划列

    def get_internal_force(self):
        """得到每一根杆件的杆端位移与杆端内力"""
        temp = []
        for i in range(len(self.in_result)):
            if self.in_result[i]:
                temp.append(i)  # temp里面存的是划行划列之前的坐标
        for element in self.elements:
            for j in range(6):
                if element.freedom[j] in temp:
                    element.delta[j] = self.result[temp.index(int(element.freedom[j]))]  # 对号入座得到整体坐标下的杆端位移
            # 局部坐标下的杆端内力
            element.force = element.T @ element.k_e
            element.force = element.force @ element.delta.T - element.T @ element.load

    def resolve(self):
        """求解行列式"""
        self.result = np.linalg.solve(self.K, self.load)


def is_splice(node, element_set: set[Element]):
    """判断该节点是否是全铰接点"""
    for element in element_set:
        if node not in element.nodes:
            raise ValueError("element not include node")
        if not element.link_way[element.nodes.index(node)]:
            return False
    return True
