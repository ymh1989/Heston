## Heston ##

### Introduction ###
The Black-Scholes model[1] concerns with the option pricing problems and has achieved great success, especially in stock option pricing. It derived a parabolic second order PDE for the value u(s,t) of an option on stocks. However, it is a well-known fact that in actual markets the Black-Scholes assumptions are violated. The most apparent violation is that the volatility implied from traded options, the implied volatility, is not constant but exhibits a strike dependency and a term structure. This strike dependency is usually referred to as the implied volatility smile. The stochastic properties of interest rate r and volatility σ often lead to more complicated approaches because analytic solutions are usually absent. A stochastic volatility model that has been particularly successful at explaining the implied volatility. It was first introduced to traditional model in the 1980’s when Hull and White[2] and others generalized the Black-Scholes model. Heston model[3] was later presented in 1993 which offered an analytic formula in semiclosed-form for the price of a vanilla option. A virtu of the Heston model is that, contrary to e.g. local volatility which is important when pricing forward skew dependent claims. It make the Heston model a prominent candidate for valuing and
hedging exotic options. Moreover, a closed-form solution exists and the calculations for Greeks are more straightforward. Heston model is a nice benchmark in testing numerical schemes dealing with parabolic partial differential equations(PDE). If there are no exists exact solution of PDE, then we can approximation by numerical computation.

In this repo, I have implemented numerical solution of Heston model[3] by using finite difference method. Especially, I would like to show how boundary condition affects the solution. I applied linear boundary condition and hybrid boundary condition is used.

###About Heston###
- `PDF` file which is described about Heston is attached in `code` folder.

### Implementation ###
- `MATLAB` codes and figures are uploaded.
- Operator spliting method(OSM) is used. 

### Future works ###
- In general, Heston model in finite difference method has to be used PDE boundary conditions at far-field area. As you can see in figure, linear boundary condtion and proposeed method are not accurate enough. Therefore, boundary condition to solmore accurate and fast is needed.


###Reference###

\[1\] Black, Fischer, and Myron Scholes. "The pricing of options and corporate liabilities." The journal of political economy (1973): 637-654.

\[2\] Hull, John, and Alan White. "The pricing of options on assets with stochastic volatilities." The journal of finance 42.2 (1987): 281-300.

\[3\] Heston, Steven L. "A closed-form solution for options with stochastic volatility with applications to bond and currency options." Review of financial studies 6.2 (1993): 327-343.