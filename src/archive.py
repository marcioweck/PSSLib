import numpy as np

from deap import tools
from operator import eq

class Archive(tools.HallOfFame):

	def __init__(self, maxsize, similar=eq):
		tools.HallOfFame.__init__(self, maxsize, similar)
	
	def update(self, population):
		"""Update the hall of fame with the *population* by replacing the
		worst individuals in it by the best individuals present in
		*population* (if they are better). The size of the hall of fame is
		kept constant.
		
		:param population: A list of individual with a fitness attribute to
                           update the hall of fame with.
		"""
		
		if len(self) == 0 and self.maxsize !=0:
			# Working on an empty hall of fame is problematic for the
			# "for else"
			self.insert(population[0])

		for ind in population:
			if ind.fitness > self[-1].fitness or len(self) < self.maxsize:
				has_twin = -1
				for i, hofer in enumerate(self):
					# Loop through the hall of fame to check for any
					# similar individual
					if self.similar(ind, hofer):
						if ind.fitness > hofer.fitness:
							has_twin = i
						break
				else:
					# The individual is unique and strictly better than
					# the worst
					if len(self) >= self.maxsize:
						self.remove(-1)
					self.insert(ind)

				if has_twin > -1:
					self.remove(has_twin)
					self.insert(ind)