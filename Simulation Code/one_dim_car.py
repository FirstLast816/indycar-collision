

import math
import matplotlib.pyplot as plt
import random


class Particle:

    def __init__(self, pos, mass, vel, order):
        self.order = order
        self.vel = vel
        self.pos = pos
        self.mass = mass
        self.traject = [[self.pos[0], self.pos[1]]]
        self.adj = []
        self.forces = []
        self.total_force = []
        self.accel = 0
        self.breaks = []

    def add_adj(self, adj):
        self.adj.append(adj)

    def get_order(self):
        return self.order

    def get_vel(self):
        return self.vel

    def get_pos(self):
        return self.pos

    def get_mass(self):
        return self.mass

    def get_traject(self):
        return self.traject

    def get_adj(self):
        return self.adj

    def get_force(self):
        return self.forces

    def get_total_force(self):
        return self.total_force

    def get_breaks(self):
        return self.breaks

    def get_adj_orders(self):
        adj_orders = [adj.get_order() for adj in self.adj]
        return adj_orders

    def update_force(self):
        self.forces = []
        for adj in self.adj:
            # Compute the magnitude of the force
            dist = math.sqrt((self.pos[0] - adj.get_pos()[0]) ** 2
                              + (self.pos[1] - adj.get_pos()[1]) ** 2)
            strain = dist / partition_length - 1
            stress = get_stress(strain)
            force_mag = (stress * 0.01 * 0.01) / (10 ** 6) # can be negative for compression
            # Compute the force vector by determining the direction
            point = [adj.get_pos()[0] - self.pos[0], adj.get_pos()[1] - self.pos[1]]
            cos_theta = point[0] / dist
            sin_theta = point[1] / dist
            force = [force_mag * cos_theta, force_mag * sin_theta]
            self.forces.append(force)
        self.total_force = [sum([self.forces[i][0] for i in range(len(self.forces))]),
                            sum([self.forces[i][1] for i in range(len(self.forces))])]

    def update_accel(self):
        if self.total_force != [0, 0]:
            force_mag = math.sqrt(self.total_force[0] ** 2 + self.total_force[1] ** 2)
            accel_mag = force_mag / mass
            cos_theta = self.total_force[0] / force_mag
            sin_theta = self.total_force[1] / force_mag
            self.accel = [accel_mag * cos_theta, accel_mag * sin_theta]
    
    def update_vel(self):
        if self.pos[0] >= wall:
            self.vel = [-1 * abs(self.vel[0]), self.vel[1]]
        self.vel[0] += self.accel[0] * time_increment
        self.vel[1] += self.accel[1] * time_increment
    
    def update_pos(self, time_incr):
        self.pos[0] += time_incr * self.vel[0]
        self.pos[1] += time_incr * self.vel[1]
        self.traject.append([self.pos[0], self.pos[1]])

    def create_cracks(self):
        for adj in self.adj:
            dist = math.sqrt((self.pos[0] - adj.get_pos()[0]) ** 2
                              + (self.pos[1] - adj.get_pos()[1]) ** 2)
            propn = dist / partition_length
            if propn < 0.929 or propn > 1.02:
                self.breaks.append(adj)
                self.adj.remove(adj)

    
class Stick:

    def __init__(self, num_parts, length, mass, pos_init, angle_init, speed_init):
        # Set properties of the stick
        self.length = length
        self.mass = mass
        vel_init = [speed_init * math.cos(angle_init), speed_init * math.sin(angle_init)]
        # Create list of particles that make up the stick
        self.particles = []
        partition_length = self.length / num_parts
        centroid_dists = []
        for i in range(num_parts):
            dist = partition_length / 2 + i * partition_length
            centroid_dists.append(dist)
        centroid_posits = []
        for i in range(num_parts):
            pos = [-1 * centroid_dists[i] * math.cos(angle_init),
                   -1 * centroid_dists[i] * math.sin(angle_init)]
            centroid_posits.append(pos)
        for i in range(num_parts):
            particle_mass = self.mass / num_parts
            order = i
            self.particles.append(Particle(centroid_posits[i], particle_mass, vel_init, order))
        self.positions = [particle.get_pos() for particle in self.particles]
        # Create empty list of cracks in the stick
        self.cracks = []
        # Create empty list of fragments in the stick
        self.fragments = []
        # Define adjacent particles for each particle in the stick
        self.particles[0].add_adj(self.particles[1])
        self.particles[num_parts - 1].add_adj(self.particles[num_parts - 2])
        for i in range(num_parts):
            if i != 0 and i != num_parts - 1:
                self.particles[i].add_adj(self.particles[i - 1])
                self.particles[i].add_adj(self.particles[i + 1])
        
    def get_length(self):
        return self.length
    
    def get_particles(self):
        return self.particles

    def get_positions(self):
        return self.positions

    def get_cracks(self):
        return self.cracks

    def get_fragments(self):
        return self.fragments

    def update_positions(self):
        self.positions = [particle.get_pos() for particle in self.particles]

    def compile_cracks(self):
        for i in range(len(self.particles) - 1):
            if self.particles[i].get_breaks().count(self.particles[i + 1]) == 1:
                self.cracks.append([self.particles[i].get_order(),
                                    self.particles[i + 1].get_order()])

    def compile_fragments(self):
        self.fragments.append(list(range(self.cracks[0][0] + 1)))
        for i in range(len(self.cracks)):
            if i == len(self.cracks) - 1:
                self.fragments.append(list(range(self.cracks[i][1], len(self.particles))))
            else:
                self.fragments.append(list(range(self.cracks[i][1], self.cracks[i + 1][0] + 1)))

    def __str__(self):
        print_str = ''
        for i in range(len(self.particles)):
            string = ('Particle ' + str(i) + '\n' + 'Mass: '
                      + str(self.particles[i].get_mass()) + ' kg.' + '\n' + 'Position: '
                      + str(self.particles[i].get_pos()) + ' m.' + '\n' + 'Velocity: '
                      + str(self.particles[i].get_vel()) + ' m./s.' + '\n'
                      + 'Adjacent Particles: ' + str(self.particles[i].get_adj_orders())
                      + '\n' + 'Forces: ' + str(self.particles[i].get_force()) + ' N.' + '\n'
                      + 'Total Force: ' + str(self.particles[i].get_total_force()) + ' N.')
            if i == len(self.particles) - 1:
                print_str += string
            else:
                print_str = print_str + string + '\n'
        return print_str


def get_stress(strain):
    # Compute stress under compression
    if strain < 0:
        compress = 1 - strain
        # Stress is negative as the particles would be pushing against each other
        stress = (907986.8024 * compress ** 3 - 89109.2718 * compress ** 2
                  + 6037.7371 * compress + 0.3962) * (-1)
        return stress
    # Compute stress under tensile forces
    elif strain < 0.0010781:
        stress = (200/0.0031) * strain
        return stress
    else:
        stress = 278.7055 / (1 + 8.6358 * math.exp(-978.5576 * strain))
        return stress


def random_color():
    digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b',
              'c', 'd', 'e', 'f']
    rand_color = '#'
    for i in range(6):
        rand_color += random.choice(digits)
    return rand_color


if __name__ == '__main__':
    # Set the parameters for the simulation
    speed_init = 70 # Average speed of an Indycar in m/s
    length = 5.1 # Length of an Indycar in m
    mass = 1.377 # Mass of aluminum beam in kg
    pos_init = [0, 0] # position in m
    partitions = 1000
    angle_init = math.pi / 3
    partition_length = length / partitions
    stick = Stick(partitions, length, mass, pos_init, angle_init, speed_init)
    # print(stick)
    time_increment = .001 # Time between position updates in s
    time = .08 # Time of simulation in s
    iterations = int(time / time_increment)
    wall = 0 # Horizontal distance from origin to wall in m


    init_positions = stick.get_positions()
    initial_x_positions = [position[0] for position in init_positions]
    initial_y_positions = [position[1] for position in init_positions]
    plt.plot(initial_x_positions, initial_y_positions, 'ko-')
    # Run the simulation
    for iteration in range(iterations):
        for particle in stick.get_particles():
            particle.update_pos(time_increment)
            particle.update_force()
            particle.update_accel()
            particle.update_vel()
        stick.update_positions()
        for particle in stick.get_particles():
            particle.create_cracks()
        # Print out data for each particle for each iteration in the simulation
        # print('\n' + "Iteration " + str(iteration + 1))
        # print(stick)
    stick.compile_cracks()
    stick.compile_fragments()
    fragments = stick.get_fragments()
    print(len(fragments))

    deg_angle = round(angle_init * 180 / math.pi)

    for fragment in fragments:
        final_positions = stick.get_positions()[fragment[0]: fragment[len(fragment) - 1] + 1]
        final_x_positions = [position[0] for position in final_positions]
        final_y_positions = [position[1] for position in final_positions]
        color = random_color()
        plt.plot(final_x_positions, final_y_positions, color=color, linestyle='-', marker='o')
    
    plt.xlabel('Distance (m.)')
    plt.ylabel('Distance (m.)')
    plt.title(str(partitions) + '-Partition Beam Crash (' + str(deg_angle) + ' Degrees)')
    # Must x out to see the next graph
    plt.show()

    partition_length = stick.get_length() / partitions
    fragment_sizes = [(len(fragment) - 1) * partition_length for fragment in fragments]
    plt.xlabel('Fragment Length (m.)')
    plt.ylabel('Frequency')
    plt.title(str(partitions) + '-Partition Beam Fragment Sizes (' + str(deg_angle) + ' Degrees)')
    plt.hist(fragment_sizes)
    # Must x out to end the program
    plt.show()
