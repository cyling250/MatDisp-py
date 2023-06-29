from Element import *


class Structure:
    def __init__(self):
        super(Structure, self).__init__()
        self.nodes = []  # 节点数组
        self.elements = []  # 单元数组
        self.size_of_K = 0  # 整体刚度矩阵维度
        self.in_result = []  # 位移定位矩阵，协助内力计算
        self.K = []  # 整体刚度矩阵
        self.load = []  # 荷载列阵
        self.result = []  # 结果矩阵

    def link(self, n1, n2, link_way1=False, link_way2=False, E=sym.Symbol("E"), I=sym.Symbol("I"),
             A=sym.Symbol("A")):
        # 杆件生成函数
        self.elements.append(Element(n1 - 1, n2 - 1, link_way1, link_way2, E, I, A))  # 生成杆件并加入杆件列表
        self.nodes[n1 - 1].elements.append([len(self.elements) - 1, link_way1])  # 将链接信息加入节点1
        self.nodes[n2 - 1].elements.append([len(self.elements) - 1, link_way2])  # 将链接信息加入节点2
        # 杆件生成完毕

    def get_L(self):
        # 计算每个单元的单元长度
        for i in range(len(self.elements)):
            L = (self.nodes[self.elements[i].nodes[0][0]].loca[0] - self.nodes[self.elements[i].nodes[1][0]].loca[
                0]) ** 2 + (
                        self.nodes[self.elements[i].nodes[0][0]].loca[1] -
                        self.nodes[self.elements[i].nodes[1][0]].loca[1]) ** 2
            self.elements[i].L = L ** 0.5

    def get_ke(self):
        # 为每个杆件生成单元刚度矩阵
        for i in range(len(self.elements)):
            E = self.elements[i].E
            L = self.elements[i].L
            A = self.elements[i].A
            I = self.elements[i].I
            temp = [
                [E * A / L, 0, 0, -E * A / L, 0, 0],
                [0, 12 * E * I / pow(L, 3), 6 * E * I / pow(L, 2), 0, -12 * E * I / pow(L, 3), 6 * E * I / pow(L, 2)],
                [0, 6 * E * I / pow(L, 2), 4 * E * I / L, 0, -6 * E * I / pow(L, 2), 2 * E * I / L],
                [-E * A / L, 0, 0, E * A / L, 0, 0],
                [0, -12 * E * I / pow(L, 3), -6 * E * I / pow(L, 2), 0, 12 * E * I / pow(L, 3), -6 * E * I / pow(L, 2)],
                [0, 6 * E * I / pow(L, 2), 2 * E * I / L, 0, -6 * E * I / pow(L, 2), 4 * E * I / L]]
            if self.elements[i].nodes[0][1] and self.elements[i].nodes[1][1]:
                # 说明这是一个轴力杆,对该轴力杆的弯矩进行自行约束
                temp = [
                    [E * A / L, 0, 0, -E * A / L, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [-E * A / L, 0, 0, E * A / L, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0]]
            temp = np.array(temp)  # 转化为数组
            # print(temp)
            cosa = (self.nodes[self.elements[i].nodes[1][0]].loca[0] - self.nodes[self.elements[i].nodes[0][0]].loca[
                0]) / L
            sina = (self.nodes[self.elements[i].nodes[1][0]].loca[1] - self.nodes[self.elements[i].nodes[0][0]].loca[
                1]) / L
            T = [
                [cosa, sina, 0, 0, 0, 0],
                [-sina, cosa, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, cosa, sina, 0],
                [0, 0, 0, -sina, cosa, 0],
                [0, 0, 0, 0, 0, 1]]
            T = np.array(T)
            self.elements[i].T = T
            Tt = np.transpose(T)
            self.elements[i].k_e = np.matmul(np.matmul(Tt, temp), T)  # 返回T转置*k*T
            # print("第{}个单元整体坐标下的单元刚度矩阵为:".format(i), self.elements[i].K_e)

    def get_freedom(self):
        # 为每个杆件生成定位向量，并编号
        freedom_num = 0  # 初始化自由度编号
        for i in self.nodes:
            link_temp = 0
            if i.is_splice():
                # 如果i是全铰接点
                for j in i.elements:
                    # 写入定位向量
                    if self.elements[j[0]].nodes[0][0] == self.nodes.index(i):  # 这个该单元的第一个节点
                        # print("单元编号",j[0])
                        self.elements[j[0]].freedom[0] = freedom_num
                        self.elements[j[0]].freedom[1] = freedom_num + 1
                        self.elements[j[0]].freedom[2] = freedom_num + 2 + link_temp
                    else:  # 这是该单元的第二个节点
                        self.elements[j[0]].freedom[3] = freedom_num
                        self.elements[j[0]].freedom[4] = freedom_num + 1
                        self.elements[j[0]].freedom[5] = freedom_num + 2 + link_temp
                    link_temp += 1
            else:
                # 如果i不是全铰接点
                add = True  # 用来判断是不是第一个结点，第一个结点不管是铰接点还是刚节点都不用+link_temp
                for j in i.elements:
                    # 写入定位向量
                    if self.elements[j[0]].nodes[0][0] == self.nodes.index(i):  # 这是该单元的第一个节点
                        if j[1]:  # 如果该节点与该单元的链接方式为铰接
                            self.elements[j[0]].freedom[0] = freedom_num
                            self.elements[j[0]].freedom[1] = freedom_num + 1
                            self.elements[j[0]].freedom[2] = freedom_num + 2 + link_temp
                            link_temp += 1
                        else:  # 如果链接方式是刚接
                            self.elements[j[0]].freedom[0] = freedom_num
                            self.elements[j[0]].freedom[1] = freedom_num + 1
                            self.elements[j[0]].freedom[2] = freedom_num + 2
                            if add:
                                # 如果是第一个结点，则需要给linktemp+1
                                link_temp += 1
                                add = False
                    else:  # 这是该单元的第二个节点
                        if j[1]:  # 如果该节点与该单元的链接方式为铰接
                            self.elements[j[0]].freedom[3] = freedom_num
                            self.elements[j[0]].freedom[4] = freedom_num + 1
                            self.elements[j[0]].freedom[5] = freedom_num + 2 + link_temp
                            link_temp += 1
                        else:  # 如果链接方式是刚接
                            self.elements[j[0]].freedom[3] = freedom_num
                            self.elements[j[0]].freedom[4] = freedom_num + 1
                            self.elements[j[0]].freedom[5] = freedom_num + 2
                            if add:
                                # 如果是第一个结点，则需要给linktemp+1
                                link_temp += 1
                                add = False
            # 定位向量书写完成，目前是没有进行去除的
            for k in range(freedom_num, freedom_num + 3 + link_temp):
                self.nodes[self.nodes.index(i)].freedom.append(k)
            freedom_num = freedom_num + 2 + link_temp

        return freedom_num

    def get_K(self):
        # 整体刚度矩阵
        self.K = [[0 for i in range(self.size_of_K)] for j in range(self.size_of_K)]
        for k in range(len(self.elements)):
            Lambda = self.elements[k].freedom  # 获取定位向量
            for i in range(6):
                for j in range(6):
                    self.K[Lambda[i]][Lambda[j]] += self.elements[k].k_e[i][j]  # 对号入座

    def set_load(self, load):
        # Load是一个二维数组，第一维度有两位[0,0],第一个表示荷载类型，True表示节点，False表示单元，第二个存储编号
        # 第二位存储着不同情况下的荷载参数形式
        for i in load:
            if i[0]:
                # 节点荷载，直接加上去就行
                self.nodes[i[1] - 1].set_load(i[2])
            else:
                # 杆件荷载，需要不同形式调用
                self.elements[i[1] - 1].set_load(i[2])

    def get_load(self):
        self.load = [0 for i in range(self.size_of_K)]
        # 先处理节点荷载
        for i in range(len(self.nodes)):
            Lambda = self.nodes[i].freedom[0:3]  # 获取定位向量
            for j in range(3):
                self.load[Lambda[j]] += self.nodes[i].load[j]  # 结点不需要反号
        # 再处理等效荷载
        for i in range(len(self.elements)):
            Lambda = self.elements[i].freedom  # 获取定位向量
            L = self.elements[i].L
            # print("局部坐标系下的单元荷载列阵:", self.elements[i].Load)
            self.elements[i].load = -np.matmul(np.transpose(self.elements[i].T),
                                               self.elements[i].load)  # 坐标转换，局部坐标转整体坐标
            # print()
            # print("整体坐标系下的单元荷载列阵:", self.elements[i].Load)
            for j in range(6):
                # print("当前正在填充第{}位荷载".format(Lambda[j]))
                self.load[Lambda[j]] += self.elements[i].load[j]

    def change_freedom(self, delete_temp=[]):
        # 对整体刚度矩阵和荷载列阵进行划行划列
        self.in_result = np.zeros(self.size_of_K)  # 初始化位移定位矩阵
        for i in range(len(self.elements)):
            for j in range(6):
                if self.elements[i].restrain[j] < 0:  # 该自由度被约束住了
                    self.in_result[self.elements[i].freedom[j]] = -1  # 位移约束赋值-1
                    delete_temp.append(self.elements[i].freedom[j])
            if self.elements[i].nodes[0][1] and self.elements[i].nodes[1][1]:
                # 说明这是一个轴力杆,对该轴力杆的弯矩进行自行约束
                self.in_result[self.elements[i].freedom[2]] = -1
                self.in_result[self.elements[i].freedom[5]] = -1
                delete_temp.append(self.elements[i].freedom[2])
                delete_temp.append(self.elements[i].freedom[5])
        if len(self.K) > 0:
            self.K = np.delete(self.K, delete_temp, axis=0)  # 约束划行
            self.K = np.delete(self.K, delete_temp, axis=1)  # 约束划列
        if len(self.load) > 0:
            self.load = np.delete(self.load, delete_temp, axis=0)  # 约束划列

    def add_restrain(self, nature, i, restrain):
        if nature:
            # 如果是节点,nature==True表示节点
            self.nodes[i - 1].add_restrain(restrain)
        else:  # 如果是杆件
            self.elements[i - 1].add_restrain(restrain)

    def get_internal_force(self):
        # 得到每一根杆件的杆端位移与杆端内力
        # 先处理位移定位矩阵
        temp = []
        for i in range(len(self.in_result)):
            if self.in_result[i] >= 0:
                temp.append(i)  # temp里面存的是划行划列之前的坐标
        for i in range(len(self.elements)):
            for j in range(6):
                if self.elements[i].freedom[j] in temp:
                    self.elements[i].delta[j] = self.result[temp.index(self.elements[i].freedom[j])]  # 对号入座得到整体坐标下的杆端位移
            # print("第{}根杆件的杆端位移为".format(i), self.elements[i].Delta)
            # 局部坐标下的杆端内力
            self.elements[i].force = np.matmul(self.elements[i].T, self.elements[i].k_e)
            self.elements[i].force = np.matmul(self.elements[i].force,
                                               np.transpose(self.elements[i].delta)) + np.matmul(
                -self.elements[i].T, self.elements[i].load)
        return
