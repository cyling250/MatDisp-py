class Node:
    def __init__(self, x=0, y=0):
        super(Node, self).__init__()
        # 节点编号就采用数组下标
        self.loca = [x, y]  # 节点坐标初始化0，0
        self.elements = []  # 节点链接的单元[单元编号，链接方式],链接方式为False是刚接，True是铰接
        self.freedom = []  # 节点自由度编号，一般有三个元素，半铰接点有特殊处理函数
        self.load = [0, 0, 0]  # 荷载列阵,节点上的荷载只考虑前三个，节点平动和节点转动

    def is_splice(self):
        # 判断节点是否是全铰接点
        for i in self.elements:
            if not i[1]:
                # 如果有一个刚接，则不是全铰接点，处理方式不需要特殊处理
                return False  # 非全铰接点
        return True  # 全铰接点

    def add_restrain(self, restrain):
        # 节点约束条件
        for i in range(3):
            if restrain[i] == -1:
                self.freedom[i] = -1
        return

    def set_load(self, load):
        # 荷载列阵
        for i in range(3):
            self.load[i] += load[i]
        return
