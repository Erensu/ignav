# ignav
integrated navigation for ins and gnss

# 简介
IGNAV基于RTKLIB开发的INS/GNSS组合导航算法库，采用C语言编写；IGNAV适用于车载场景，目前正进行INS/GNSS融合视觉信息、激光雷达的算法编写，后续会逐步更新。
IGNAV主要功能包括：
  1. 松耦合算法；
  2. SPP、PPP、DGPS、RTK紧耦合算法；
  3. 里程计辅助；
  4. 磁力计辅助；
  5. NHC、ZUPT、ZARU等运动约束；
  6. Doppler辅助INS/GNSS；
  7. 双天线航向辅助；
  8. 静对准、动对准初始化；
  9. INS正向和反向机械编排；
  10. INS/GNSS正向和反向组合滤波；
  11. 初步支持视觉信息辅助定位定姿；
  12. RTS/前后向滤波平滑；
  13. 车载轨迹动态显示等。

# 配置
IGNAV在解算数据前，需要编辑相应的配置文件，配置文件的编写可以参考RTKLIB的说明文档，配置文件在/example/conf目录。

# 解算
默认IGNAV解算结果放在bin文件navapp所在目录。

# 编译
IGNAV目前仅支持在Linux环境下编译，具体过程与编译CMAKE工程一样。

# 运行
在目录/example/demo下有相应视频Demo；
plot运行需要安装Qt．
运行命令：navapp -o ../example/conf/navlib.conf -m 52716

# 参考
[1] http://www.rtklib.com

[2] P.D.Groves,Principles of GNSS,Intertial,and Multisensor Integrated Navigation System

[3] 武汉迈普导航科技有限公司，http://www.whmpst.com/cn/

# 联系
地址：武汉大学，卫星导航与定位技术研究中心，苏景岚

QQ  ：1971129844

邮箱：sujinglan@whu.edu.cn

