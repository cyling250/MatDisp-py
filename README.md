# 矩阵位移法程序说明文档
&emsp;&emsp;本程序的编制目的是为了对华中科技大学出版社出版的《结构力学》教材（下称教材）第九章矩阵位移法进行程序实现。
## 程序流程图
&emsp;&emsp;本程序采用PYTHON语言进行程序设计，有效利用PYTHON中提供的类与对象方法，创建了结点类Node、杆件单元类Element和整体结构类Structure。这三个类的创建使得本程序可以方便求解如梁、刚架、桁架、组合结构等平面杆系结构的内力和位移，并支持符号运算。因此，本程序具有较好的通用性。
<div align=center>程序流程图</div>
</br>
<div align=center><img src="flow_chart.png"></div>
&emsp;&emsp;上图为程序计算流程，具体步骤如下：</br>
&emsp;&emsp;1）定义结点、杆件，并生成结构。连接杆件完成后，程序可以自行组合成整体结构。</br>
&emsp;&emsp;2）对每根杆件生成局部荷载列阵和局部刚度矩阵。</br>
&emsp;&emsp;3）使用前处理法进行自由度编号，生成每根杆件的定位向量。</br>
&emsp;&emsp;4）依据定位向量对号入座，生成整体刚度矩阵和整体荷载列阵。</br>
&emsp;&emsp;5）读取边界约束条件，对整体刚度矩阵和整体荷载列阵进行划行划列。</br>
&emsp;&emsp;6）联立求解基本方程，得到杆端位移和杆端内力。</br>

## 数据定义

### Node类
&emsp;&emsp;构造函数Node(x,y)
#### 1.数据成员
&emsp;&emsp;loca：结点的坐标，用一个定长为2的列表表示，列表第一个元素为结点x坐标，第二个元素为结点y坐标。</br>
&emsp;&emsp;elements：结点所链接的单元，用变长列表表示，元素存储格式为[单元编号,链接方式]。全局定义链接方式False表示刚接，True表示铰接。</br>
&emsp;&emsp;freedom：结点自由度编号，用边长列表表示，一般长度为3，前三个自由度分别表示x方向平动、y方向平动和转动。当结点有铰接时，结点的自由度可能会大于3.</br>
&emsp;&emsp;load：结点荷载列阵，用定长为3的列表表示，依次表示x方向荷载，y方向荷载和弯矩荷载。</br>

#### 2.函数成员
&emsp;&emsp;bool is_splice(None)：判断当前结点是否是全铰接点，如果是全铰接点返回True，否则返回False。</br>
&emsp;&emsp;void add_restrain(restrain):依据传入的restrain数据，对结点自由度进行赋值，被约束住赋值为-1。</br>
&emsp;&emsp;void set_load(load):依据传入的load数据，对结点荷载进行赋值。</br>
### Element类
&emsp;&emsp;构造函数Element(node1,node2,link_way1,link_way2,E,I,A)
#### 1.数据成员
&emsp;&emsp;nodes：定长2×2列表，存放单元所包含的两个结点编号及其链接方式，元素结构为[node,link_way]。</br>
&emsp;&emsp;freedom：单元定位向量及其自由度，定长为6的列表，元素不为-1时，描述为单元的定位坐标，元素为-1时，表示该自由度被约束。</br>
&emsp;&emsp;E：单元弹性模量。</br>
&emsp;&emsp;I：单元截面惯性矩。</br>
&emsp;&emsp;L：单元长度。</br>
&emsp;&emsp;A：单元截面面积。</br>
&emsp;&emsp;k_e：整体坐标系下的单元刚度矩阵。</br>
&emsp;&emsp;T：单元的坐标变换矩阵。</br>
&emsp;&emsp;restrain:单元的约束条件，定长为6的列表，元素为-1时表示该自由度被约束。</br>
&emsp;&emsp;load：局部坐标下荷载列阵。</br>
&emsp;&emsp;delta：杆端位移，为后处理矩阵。</br>
&emsp;&emsp;force：杆端内力，为后处理矩阵。</br>
#### 2.函数成员
&emsp;&emsp;void add_restrain(restrain):依据传入的restrain对单元施加约束条件。</br>
&emsp;&emsp;void set_load(load):依据传入的load对单元施加荷载。输入参数load为一个不确定形状的列表：当load[0]为0时，表示传入了杆上集中荷载，此时load[1]表示集中荷载大小，load[2]表示集中荷载到杆件A端的距离；当load[0]为1时，表示杆端弯矩，此时load[1]表示作用在杆件A端的弯矩，load[2]表示作用在杆件B端的弯矩;当load[0]为2时，表示均布荷载，此时load[1]表示均布荷载的集度；当load[0]为3时，表示三角形均布荷载，此时load[1]表示三角形均布荷载最大值。荷载的定义方式可以自行定义，只要最终将荷载处理到单元的数据成员load上即可。
### structure类
#### 1.数据成员
&emsp;&emsp;nodes：存放真实结点的数组，一个元素为一个结点对象，该数组的下标为对应结点的编号。</br>
&emsp;&emsp;elements：存放真实单元的数组，一个元素为一个单元对象，该数组的下标为对应单元的编号。</br>
&emsp;&emsp;size_of_K：整体刚度矩阵的维度。</br>
&emsp;&emsp;in_result：位移的定位矩阵，用来后处理协助内力计算。</br>
&emsp;&emsp;K：整体刚度矩阵。</br>
&emsp;&emsp;load：整体荷载列阵。</br>
&emsp;&emsp;result：位移矩阵，为求解整体刚度方程的结果。</br>
#### 2.函数成员
&emsp;&emsp;void link(n1,n2,link_way1,link_way2,E,I,A)：杆件单元生成函数，用来生成每一根杆件。</br>
&emsp;&emsp;void get_L(None)：内部调用函数，用来生成每一根杆件的长度。</br>
&emsp;&emsp;void get_ke(None)：用来生成每一根杆件的局部刚度矩阵。</br>
&emsp;&emsp;void get_freedom(None)：用来生成每一根杆件的定位向量。</br>
&emsp;&emsp;void get_K(None)：用来组合整体刚度矩阵。</br>
&emsp;&emsp;void set_load(load)：用来将荷载作用到结构上，但是该函数不会生成荷载列阵。</br>
&emsp;&emsp;void get_load(None)：用来生成结构荷载列阵。</br>
&emsp;&emsp;void change_freedom(delete_temp)：用来实现约束条件的施加，划行划列。</br>
&emsp;&emsp;void add_restrain(nature,i,restrain)：用来施加约束条件，nature表示约束条件的类型，nature为True表示为结点约束，nature为False表示为杆端约束。</br>
&emsp;&emsp;void get_internal_force(None)：后处理函数，用来生成每一根杆件的杆端内力。</br>
&emsp;&emsp;本程序使用教材例9-15进行了验证，但难免有其他方面的错误与不足，恳请批评指正！

&emsp;&emsp;此外，还提供了c++实现矩阵位移法的程序，c++相比于python具有更高的运行速度和清晰度，但代码结构整体比python复杂，以供对比学习。详见</br>
&emsp;&emsp;https://github.com/cyling250/MSDinSM-cpp</br>