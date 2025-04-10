# Quantum Mechanical Scattering Problems / 量子力学散射问题

## Introduction / 简介

This repository contains a collection of programs for calculating differential and integral cross sections for quantum mechanical scattering problems from reactance or transition matrices. The main components include:

本代码库包含一系列用于从反应矩阵或跃迁矩阵计算量子力学散射问题的微分和积分截面的程序。主要组件包括：

- **SECTIONS_FOR_QUANTUM_MECHANICAL_SCATTERING_PROBLEMS_FROM_REACTANCE.f**: The main Fortran program for calculating differential and integral cross sections from reactance or transition matrices.
  
  用于从反应矩阵或跃迁矩阵计算微分和积分截面的主要Fortran程序。

- **Frescox**: A coupled-channels calculation code for nuclear reactions (version 7.2.2).
  
  用于核反应的耦合通道计算代码（7.2.2版本）。

- **WavePacket**: A MATLAB package (version 5.3.0) for quantum-mechanical wavepacket dynamics simulations.
  
  用于量子力学波包动力学模拟的MATLAB软件包（5.3.0版本）。

- **Ion-Atom-Wave**: Program for calculation of single ionization cross sections in ion-atom collisions.
  
  用于计算离子-原子碰撞中单电离截面的程序。

- **Machine-Learning-Molecular-Dynamics**: Machine learning approaches for molecular dynamics simulations.
  
  用于分子动力学模拟的机器学习方法。

## Main Program / 主程序

The main program `SECTIONS_FOR_QUANTUM_MECHANICAL_SCATTERING_PROBLEMS_FROM_REACTANCE.f` calculates the differential and integral cross sections using formulas from J. Blatt and L. Biedenharn (Rev. Mod. Phys. 24, 258, 1952) with corrections by R. Huby (Proc. Phys. Soc. London, 67A, 1103, 1954).

主程序`SECTIONS_FOR_QUANTUM_MECHANICAL_SCATTERING_PROBLEMS_FROM_REACTANCE.f`使用J. Blatt和L. Biedenharn（Rev. Mod. Phys. 24, 258, 1952）的公式计算微分和积分截面，并采用了R. Huby（Proc. Phys. Soc. London, 67A, 1103, 1954）的修正。

### Features / 特性:

- Calculation of differential and integral cross sections
  
  计算微分和积分截面

- Processing of reactance (R) or transition (T) matrices
  
  处理反应矩阵(R)或跃迁矩阵(T)

- Handling of various quantum mechanical scattering systems
  
  处理各种量子力学散射系统

## Frescox / Frescox程序

Frescox is a scattering code for coupled-channels calculations in nuclear physics. This is version 7.2.2 available from [https://github.com/LLNL/Frescox](https://github.com/LLNL/Frescox).

Frescox是用于核物理耦合通道计算的散射代码。这是7.2.2版本，可从[https://github.com/LLNL/Frescox](https://github.com/LLNL/Frescox)获取。

## WavePacket / WavePacket程序

WavePacket is a program package for quantum-mechanical simulations. This repository contains version 5.3.0 released on May 30, 2017.

WavePacket是一个用于量子力学模拟的程序包。本代码库包含2017年5月30日发布的5.3.0版本。

### Installation / 安装:

Add the Sources directory to your MATLAB source path. For each simulation, write an initialization function `qm_init()`, change the working directory to where this function is located, and run:

将Sources目录添加到你的MATLAB源代码路径。对于每次模拟，编写一个初始化函数`qm_init()`，将工作目录更改为该函数所在的位置，然后运行：

```matlab
qm_setup();
qm_init();
qm_propa();
qm_cleanup();
```

## Requirements / 需求

- Fortran compiler (for the main program and Frescox)
  
  Fortran编译器（用于主程序和Frescox）

- MATLAB (for WavePacket)
  
  MATLAB（用于WavePacket）

## License / 许可证

The included software components are subject to their respective licenses:

包含的软件组件受各自许可证的约束：

- Frescox: GNU General Public License v2 or later
  
  Frescox：GNU通用公共许可证v2或更高版本

- WavePacket: GNU General Public License v2 or later
  
  WavePacket：GNU通用公共许可证v2或更高版本

- Other components may have their own licenses, please check the respective directories.
  
  其他组件可能有自己的许可证，请查看相应的目录。

## References / 参考文献

- Brandt, M.A., Truhlar, D.G., Smith, R.L. "Program for calculating differential and integral cross sections for quantum mechanical scattering problems from reactance or transition matrices." Comp. Phys. Commun. 5 (1973) 456.
  
- Blatt, J. and Biedenharn, L., Rev. Mod. Phys. 24, 258 (1952)
  
- Huby, R., Proc. Phys. Soc. (London) 67A, 1103 (1954)