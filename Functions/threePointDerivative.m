function derivative = threePointDerivative(xiMinus, yiMinus, xi, yi, xiPlus, yiPlus, pointLocation)
%From NUMERICAL METHODS FOR ENGINEERS AND SClENTISTS by Gilat and
%Subramaniam

switch pointLocation
    case 'First'
        derivative = yiMinus*((2*xiMinus - xi - xiPlus) / ((xiMinus-xi) * (xiMinus-xiPlus))) + yi*((xiMinus - xiPlus) / ((xi-xiMinus) * (xi-xiPlus))) + yiPlus*((xiMinus - xi) / ((xiPlus-xiMinus) * (xiPlus-xi)));
    case 'Interior'
        derivative = yiMinus*((xi - xiPlus) / ((xiMinus-xi) * (xiMinus-xiPlus))) + yi*((2*xi - xiMinus - xiPlus) / ((xi-xiMinus) * (xi-xiPlus))) + yiPlus*((xi - xiMinus) / ((xiPlus-xiMinus) * (xiPlus-xi)));
    case 'Last'
        derivative = yiMinus*((xiPlus - xi) / ((xiMinus-xi) * (xiMinus-xiPlus))) + yi*((xiPlus - xiMinus) / ((xi-xiMinus) * (xi-xiPlus))) + yiPlus*((2*xiPlus - xiMinus - xi) / ((xiPlus-xiMinus) * (xiPlus-xi)));
end