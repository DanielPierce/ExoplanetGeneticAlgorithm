

def Lerp(input, lowerBound, upperBound):
    return (input * lowerBound) + ((1-input) * upperBound)

def Clamp(input, lowerBound=-float('inf'), upperBound=float('inf')):
    return max(lowerBound, min(input, upperBound))