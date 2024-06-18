# 矩阵位移法程序说明文档

&emsp;&emsp;本程序的编制目的是为了对华中科技大学出版社出版的《结构力学》教材（下称教材）第九章矩阵位移法进行程序实现。

## 程序流程图

&emsp;&emsp;本程序采用PYTHON语言进行程序设计，有效利用PYTHON中提供的类与对象方法，创建了结点类Node、杆件单元类Element和整体结构类Structure。这三个类的创建使得本程序可以方便求解如梁、刚架、桁架、组合结构等平面杆系结构的内力和位移，并支持符号运算。因此，本程序具有较好的通用性。
<div align=center>程序流程图</div>
</br>
<div align=center><img src="res/flow_chart.png"></div>
&emsp;&emsp;上图为程序计算流程，具体步骤如下：</br>
&emsp;&emsp;1）定义结点、杆件，并生成结构。连接杆件完成后，程序可以自行组合成整体结构。</br>
&emsp;&emsp;2）对每根杆件生成局部荷载列阵和局部刚度矩阵。</br>
&emsp;&emsp;3）使用前处理法进行自由度编号，生成每根杆件的定位向量。</br>
&emsp;&emsp;4）依据定位向量对号入座，生成整体刚度矩阵和整体荷载列阵。</br>
&emsp;&emsp;5）读取边界约束条件，对整体刚度矩阵和整体荷载列阵进行划行划列。</br>
&emsp;&emsp;6）联立求解基本方程，得到杆端位移和杆端内力。</br>

## 验证

&emsp;&emsp;本程序使用教材例9-15进行了验证，验证程序详见[main.py](test/main.py)

## 更多

&emsp;&emsp;此外，还提供了c++实现矩阵位移法的程序，c++相比于python具有更高的运行速度和清晰度，但代码结构整体比python复杂，以供对比学习。详见</br>
&emsp;&emsp;https://github.com/cyling250/MSDinSM-cpp</br>

## [MIT-LICENSE](LICENSE)
