# LC Toolbox: its purpose

The LC Toolbox (Linear Control Toolbox) is designed to support the control engineer throughout the entire controller design procedure, i.e. starting from identification up to the design itself and the validation. This is done by providing an intuitive syntax to formulate the design in a natural way, shielding the user from the complex mathematics happening behind the scenes. The toolbox provides 3 modules: the basic systems theory module, the identification module and the controller design module.

# Installation

The LC Toolbox runs within MATLAB and requires the Robust Control Toolbox. This is a minimal setup which enables the user to use the systems theory module.

Recommended external packages:
* [YALMIP](https://yalmip.github.io/): some solvers require yalmip for the parsing of LMIs
* A decent LMI solver: this might improve the quality of a solution as well as the speed at which it is solved.
  1. [Mosek](https://www.mosek.com/): free academic license
  2. [SeDuMi](http://sedumi.ie.lehigh.edu/): free LMI solver
  3. [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html): free LMI solver
* [Optispline](https://github.com/meco-group/optispline): enables LPV controller design

# User guide (crash course)

## 1. Systems theory module

The systems theory module is the core of the toolbox and allows you to describe the plant. In order to do so, you need 3 ingredients: **models**, **systems** and **signals**.

### 1.1 Models vs Systems

The LC Toolbox makes a clear distinction between models and systems. Systems refer to a physical entity, e.g. a robot, an electric circuit, a pump, ... Each system has its own behavior which you can **try** to describe with different models. In other words, you can have multiple models for one single system.

The toolbox supports a bunch of models. 

* SSmod/DSSmod: counterparts for matlab's standard `ss` and `dss`
* ODEmod: general nonlinear dynamics of the form `dx/dt = f(x,u), y = g(x,u)`
* FRDmod: measured frequency response function
* Umod: uncertain model with nominal model and uncertainty quantification
* Gridmod: series of models, usefull for elementwise operations

To keep track of a model within the toolbox, you are encouraged to set the property `name` to the adequate value.

Systems are created using the `IOSystem(nin,nout)` call. Once a system is created, you can add as many models to it as you want using the function `add`:

```matlab
mod1 = SSmod(A,B,C,D,Ts);
mod1.name = 'linmod';
mod2 = FRDmod(freqs,resp);
mod2.name = 'measmod';

G = IOSystem(1,1);
G.add(mod1,mod2);
```

### 1.2. Signals

Once systems start interacting, you need signals to connect them. Each system has an input `in` and output `out` (which might have length > 1). We can use these to make a connection list and the `==` statement. Then, we can make a new system which represents the connected system-of-systems, as shown in the code below. P is the connected system which will have 2 inputs (S1.in(1) and S2.in(1)) and 1 output (s2.out).

```matlab
S1 = IOSystem(1,1);
S2 = IOSystem(2,1);
connections = [S1.out == S2.in(2)];
P = IOSysten(S1,S2,connections);
```

Of course, more elaborate connections can be made fairly easily. Suppose we want to make an error feedback loop. It would look as follows:

```matlab
G = IOSystem(1,1);
K = IOSystem(1,1);

r = Signal(); % we make a new signal
e = r - G.out;
connections = [K.in == e; G.in == K.out];
CL = IOSysten(G,K,connections);
```

As you can see, we needed a new Signal which we can construct via the `Signal()` constructor. Also operations on signals are supported, e.g. `2*G.out` or `G.in+G.out`. Signals can also be used to access transfer functions of the system. Suppose we want to access the transfer function going from r to K.out. This done via `CL(K.out,r)`. Note that also the internal signals remain available for inspection. No information gets lost!

### 1.3. Parameters

LCToolbox also supports Linear Parameter Varying (LPV) controller design. To this end, the SchedulingParameter class was introduced. Although it might sound complicated, it is just a parameter which you can use within your matrices. The only thing you have to do is supply a domain in which the value will remain, e.g. [-1,1] and possibly the maximal change over time, e.g. [-0.1,0.1]. Have a look at the code below: 

```matlab
range = [-1,2];
rate = [-3,2];
p = SchedulingParameter('p',range,rate);
A = [1,p^2-3;p 0];
B = [1;0];
C = [1,1];
D = 0;

mod = SSmod(A,B,C,D,p);
```

So we can use the parameter to make parameter dependent matrices and use them to make an LPV system. Note that when making a parameter dependent model (this also counts for an ODEmod!), the parameter should be added to the list of arguments.

## 2. Identification module

Coming soon! The identification module is currently being thoroughly revised and refactored to make it easier for non-experienced users. 

## 3. Controller design module

The toolbox focuses on the H-infinity/H-2 formalism which is an optimal controller design strategy. The ingredients to optimal controller design are channels, weights and norms. These are again provided by the toolbox and are readily accessible.

### 3.1. Channels

A channel is nothing more than the relation between an input and an output, for instance looking from `r` to `e`, or denoted as a transfer function: `e/r`. The toolbox uses exactly this notation to create a channel: `T = G.out/r;`. In order to keep track of transfer functions, one can also add a name to the channels using the call: `T = Channel(G.out/r,'Complementary Sensitivity')`. Channels can also be used to access the specific transfer function in a system object: `Tmod = T(CL);` will get you the model of the channel T within the closed loop. In this case, you could consider the channel to be some kind of operator.

### 3.2. Weights

Weights reflect the user's design criteria. Although any (stable!) transfer function can be used as a weight, usually only a few different weights are employed. Therefor, the Weight class is implemented. It offers these standard weights and makes them easily accessible to the user. An overview:

* Weight.DC(max): weight to constrain the maximum of a transfer function to max (dB)
* Weight.LF(crossover,order,dcgain): weight to push low frequencies down, useful to shape for instance a sensitivity function
* Weight.HF(crossover,order,dcgain): weight to push high frequencies down, useful to shape for instance a complementary sensitivity function

### 3.3. Norms

Channels and Weights are combines into a norm object. Suppose we have a channel C and a weight W, then the product of the two will result in a Norm object: `N = W*C`. This will by default be an infinity norm. If you want to use 2-norms, use the call `N = Norm(W*C,2)`.

### 3.4. Computing a controller

The important call when it comes to computing a controller is `solve`. `solve` is called on a system and takes objectives and constraints as input arguments. Suppose we want to design a controller for the standard error feedback loop CL. As objective, we have `Ws*S` and as constraints we have `Ms*S <= 1` and `Wt*T <= 1`. The only thing we have to do then is to call:

```matlab
obj = Ws*S;
constr = [Ms*S <= 1, Wt*T <= 1];
[CL,controller,info] = CL.solve(obj,constr,K);
```

If a controller exists, it will be added to K which is considered to be the optimization variable. If you are interested in some visual feedback, try `bodemag(info)` or `sigma(info)`. Note that calling T(CL)  now yields actual models as the closed loop can now be determined (it was not possible before because K was still empty!).

## 4. Simulation

You can easily simulate a solution with the toolbox. The function you need is `sim` and can be called on systems as well as on models. provide a function for the input u(t) and if used a function for the parameter p(t). Pass it on to sim and you get the simulation result!

```matlab
u = @(t) min(max(0,10*(t-1)),1); % rate limited step at time instance 1
sim(T(CL),u,10); % simulation of T using input u(t) up to t=10s.
``` 

## 5. In the pipeline

* ...