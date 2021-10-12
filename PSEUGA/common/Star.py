


class Star:
    
    
    def __init__(self, radius, mass, flux, distance):
        self.radius = radius
        self.mass = mass
        self.flux = flux
        self.distanceTo = distance

    def PrettyPrint(self):
        print(f"Radius: {self.radius} km")
        print(f"Mass: {self.mass} kg")
        print(f"Base Flux: {self.flux} e/s")
        print(f"Distance: {self.distanceTo} km")
    