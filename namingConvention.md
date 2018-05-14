# Naming convention

To ensure a similar appearance of all codes and to allow for an easier and more intuitive understanding of the purpose of identifiers (e. g., variable or function names), we agreed on the following naming convention for the source code of *FESTUNG*:

## General Rules
Two main rules should always be kept in mind: WORM and Camel Case.

### WORM
There is a very simple acronym that governs code development:

> WORM: Write Once Read Many

This means, most code parts are written once and stay like this for a long time - but are read many times to understand what is happening there. Hence, code should always be written in a way that makes it most easy for the reader to understand it, even if it means additional work for the programmer.

In clear words:
* Use **meaningful names** for functions and variables. Bytes are cheap these days, so use more than one letter to name variables and functions.
* **Structure** your code clearly. Use indentation to make it more readable and insert empty lines to form logical groups of statements. Write every statement in a new line.
* **Document** and **comment** your code. Provide a clear description for each routine and add short comments inside your algorithm to ease understanding.

If you adhere to these rules, you will not only make it easier for new people to understand the code but also help yourself when you are trying to figure out what you did there after some time.

### Camel Case
We stick to the so-called *camel case* notation. This means, words in identifiers start with an upper case letter and the rest of the word is written in lower case. Multiple words in the identifier are concatenated without any additional characters in between (no underscores!).

For **variable names and function names** we use the *lower camel case* notation, i. e., the first letter of the name is always lower case, despite it being the beginning of a word. All following words start with an upper case letter. 

Some examples:

    phi(...)
    integrateRefElemPhiPhi(...)
    dataDisc
    tEnd
    eta

For **class names** we use the *upper camel case* notation, i. e., the first letter of the name is always upper case. Since no classes are implemented, yet, this is for future code versions. 

Some examples:

    TimeStepper
    Quadrature

## Variable Names
Use the *lower camel case* notation and use meaningful names. One letter variables can be sufficient (e. g., when abbreviating the number of basis functions as ```N``` all the way through the code) but should be the exception.

### Global Variables
Global variables can be tricky and limit the general applicability of your code. When used nevertheless, they should start with the letter ```g```. For example ```gPhi2D``` are the precomputed values of the two dimensional basis functions in the quadrature points.

### Function Handles
In our numerical setting, continuous functions are rather exceptional than the rule. Continuous functions should therefore be clearly identifiable as such. We use the leading keyword ```func...``` or the suffix ```...Cont``` to make this obvious. In functions where only one function handle is passed as an argument, we use ```funcCont```.

### Coefficient Functions/Data
In *FESTUNG*, unknowns and spatially varying coefficients are commonly stored by their local (i. e. elementwise) coordinates of a globally discontinuous approximation space that uses *modal* basis functions on each element.  Modal bases do not satisfy a Lagrangian/nodal property, which is, e. g., required for visualization routines though.  *FESTUNG* provides routines for basis transformations such as [projectDataDisc2DataLagr.m](projectDataDisc2DataLagr_8m.html).  
To distinguish between modal and nodal bases, we use the suffices ```...Disc``` and ```...Lagr```, respectively. 

For example, after an L2-projection of a function *d(x,y)* into a discontinuous, modal space, we could name the array holding the coefficients ```dDisc```. In functions that expect only one coefficient argument, we may it ```dataDisc``` instead.

## Function Names
Functions always *do* something, which means we can probably describe them by a verb. This verb should always form the first part of the function name, for example as in [visualizeGrid](visualizeGrid_8m.html).

### Integration Routines
Integration routines compute blocks from integrals on the reference elements and are always named with the leading keyword `integrate...`, which is followed by the geometric entity over which we integrate - for example `integrateRefEdge...` or `integrateRefElem...`. The remainder of the function name should describe the argument of the integral. For example, the integration of the product of the two basis functions *phi* on the reference element is named `integrateRefElemPhiPhi(...)`. Derivatives of the basis function are called `Dphi`.

On edges, integrals can contain basis functions from either of the adjacent element. To distinguish those cases, we use `PhiInt` or `PhiExt` to identify the basis functions origin as the interior or exterior of the element. For example, `integrateRefEdgePhiIntPhiExt(...)` integrates over all reference edges the product of the basis function from the inside of the corresponding element and the basis function from the neighboring element.

### Assembly Functions
All assembly functions start with the verb `assemble...` and are followed by the type of the assembled array, for example `assembleMat...()` or `assembleVec...()` and the geometric entity (`...Elem...` or `...Edge...`).  The remainder should again describe the integral argument and can be built up from the following:
* `Phi`/`Dphi`: A basis function or its derivative from the reference element and a neighbouring element. Both contributions are summed up and divided by 2.
* `PhiInt`/`DphiInt`: A basis function or its derivative from the reference element.
* `PhiExt`/`DphiExt`: A basis function or its derivative from the neighbouring element.
* `FuncDisc`: A coefficient function in the discontinuous space. All integral terms are summed up (i.e., the discontinuous function value is evaluated).
* `FuncCont`: A continuous coefficient function.
* `Nu`: The integral contains the (edge) normal

