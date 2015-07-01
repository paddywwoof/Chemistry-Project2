import numpy as np

class SearchCube(object):
    """random displacement with bounds"""
    def __init__(self, xmin, xmax, interaction_manager , stepsize=0.5):
        self.xmin = xmin
        self.xmax = xmax
        self.stepsize = stepsize
        self.interaction_manager = interaction_manager

    def __call__(self, x):
        x = x.reshape((self.interaction_manager.number_signals, 3))
        """take a random step but ensure the new position is within the bounds"""

        for atom1_index,atom1 in enumerate(x):
            if min([self.calculate_magnitude(atom1-atom2) for atom2 in x if atom1 is not atom2 ]) > 50:
                atom1_interactions = self.interaction_manager.interaction_matrix[atom1_index]
                for index in atom1_interactions:
                    if index == 0:
                        continue
                    else:
                        x[atom1_index] = np.copy(x[index])
                        break




        xnew = x + np.random.uniform(-self.stepsize, self.stepsize, np.shape(x))



        return xnew
    def calculate_magnitude(self,vec):
        if vec.tolist() == [0,0,0]:
            return 0
        try:
            distance = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        except Exception as e:
            print("Vector:"+str(vec))
            raise e
        return round(distance)
