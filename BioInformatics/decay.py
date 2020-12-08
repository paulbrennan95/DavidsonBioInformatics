init_mass = int(raw_input("What is the initial mass?"))
curr_mass = int(raw_input("What is the current mass?"))
half_life = float(raw_input("What is the half life of the isotope in minutes?"))

import math

decay_rate = math.log(0.5)/half_life

decay_time = math.floor(math.log(curr_mass/float(init_mass))/decay_rate)

print "The decay rate is:", decay_rate, "\nThe decay time was:", decay_time , "minutes"