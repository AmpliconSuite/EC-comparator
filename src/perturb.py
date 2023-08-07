"""
Perturbations strategies to generate a large spectrum of different structures
which are inferred from the same structure
"""
import math
import pandas as pd
from copy import deepcopy

from utils import HEADER as h


def cartesian2polar(x, y):
	# return r and theta(degrees)
	r = (x ** 2 + y ** 2) ** .5
	theta = math.degrees(math.atan2(y, x))
	return r, theta


def polar2cartesian(r, theta):
	"""

	Arguments:
		r:
		theta: in degrees

	Return:
	"""
	x = r * math.cos(math.radians(theta))
	y = r * math.sin(math.radians(theta))
	return x, y


def translate(x, y, direction=1, step=100):
	"""
	Move coordinates to left or right
	(direction +1 or -1) using a step.

	Arguments:
		x (float):
		y (float):
		direction (int): {+1,-1}
		step (int):

	Returns:
		Coordinates (x,y)
	"""
	return x + direction * step, y + direction * step


def translate_by_radius(x, y, radius_step=10000):
	"""
	Keep teta fix and change radius
	"""
	r, theta = cartesian2polar(x, y)
	return polar2cartesian(r + radius_step, theta)


def translate_by_phi(x, y, grade=10):
	r, theta = cartesian2polar(x, y)
	return polar2cartesian(r, theta + grade)


def translate_by_radius_phi(x, y, radius_step=10000, grade=10):
	"""
	Keep teta fix and change radius
	"""
	r, theta = cartesian2polar(x, y)
	return polar2cartesian(r + radius_step, theta + grade)


def translate_by_resize(x, y, scale=2):
	return x - abs(x - y) * scale, y + abs(x - y) * scale


def extend_x(x, y, direction=1, step=100):
	"""
	Extend x coordinate in the direction
	+1 (right) or -1 (left) using the step.

	Arguments:
		x (float):
		y (float):
		direction (int): {+1,-1}
		step (int):

	Returns:
		Coordinates (x,y)
	"""
	return x + direction * step, y


def extend_y(x, y, direction=1, step=100):
	"""
	Extend y coordinate in the direction
	+1 (right) or -1 (left) using the step.

	Arguments:
		x (float):
		y (float):
		direction (int): {+1,-1}
		step (int):

	Returns:
		Coordinates (x,y)
		"""
	y_new, x_new = extend_x(y, x, direction=direction, step=step)
	return x_new, y_new


def perturb_df_translate(df_t, step):
	"""
	Translate entire data.frame
	"""
	df_new = deepcopy(df_t)
	df_new[h.START] = df_new.apply(lambda x: translate(x[h.START], x[h.END], direction=1, step=step)[0], axis=1)
	df_new[h.END] = df_new.apply(lambda x: translate(x[h.START], x[h.END], direction=1, step=step)[1], axis=1)
	return df_new


def perturb_df_scale(df_t, scale):
	"""
	Resize the cycle
	"""
	df_new = deepcopy(df_t)
	df_new[h.START] = df_new.apply(lambda x: translate_by_resize(x[h.START], x[h.END], scale=scale)[0], axis=1)
	df_new[h.END] = df_new.apply(lambda x: translate_by_resize(x[h.START], x[h.END], scale=scale)[1], axis=1)
	return df_new


def perturb_df_extend_y(df_t, step):
	"""
	Extend end position
	"""
	df_new = deepcopy(df_t)
	df_new[h.END] = df_new.apply(lambda x: extend_y(x[h.START], x[h.END], direction=1, step=step)[1], axis=1)
	return df_new


def perturb_df_extend_x(df_t, step):
	"""
	Extend start position
	"""
	df_new = deepcopy(df_t)
	df_new[h.START] = df_new.apply(lambda x: extend_x(x[h.START], x[h.END], direction=1, step=step)[0], axis=1)
	return df_new


def perturb_df_rotate(df_t, grade):
	"""
	Keep distance between the to points constant but rotate them
	"""
	df_new = deepcopy(df_t)
	df_new[h.START] = df_new.apply(lambda x: translate_by_phi(x[h.START], x[h.END], grade=grade)[0], axis=1)
	df_new[h.END] = df_new.apply(lambda x: translate_by_phi(x[h.START], x[h.END], grade=grade)[1], axis=1)
	return df_new


def perturb_df_rotate_2(df_t, step, grade):
	"""
	Keep distance between the to points constant but rotate them
	"""
	df_new = deepcopy(df_t)
	df_new[h.START] = df_new.apply(
		lambda x: translate_by_radius_phi(x[h.START], x[h.END], radius_step=step, grade=grade)[0], axis=1)
	df_new[h.END] = df_new.apply(
		lambda x: translate_by_radius_phi(x[h.START], x[h.END], radius_step=step, grade=grade)[1], axis=1)
	return df_new


def perturb_input_translate(t_file, step=100, n=10):
	"""
		Read input true and reconstruct file
		"""
	dtype = {'#chr': 'str',
			 'start': 'int',
			 'end': 'int',
			 'circ_id': 'str',
			 'estimated_cn': 'double',
			 'strand': 'str'}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)

	i = 0
	while (i < n):
		t_perturb = perturb_df_translate(t_collection, i * step)
		yield t_collection, t_perturb
		i += 1


def perturb_reconsturction_translate(t_file, r_file, step=100, n=10):
	"""
	Perturb reconstruction file by translating
	"""
	dtype = {'#chr': 'str',
			 'start': 'int',
			 'end': 'int',
			 'circ_id': 'str',
			 'estimated_cn': 'double',
			 'strand': 'str'}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)
	r_collection = pd.read_csv(r_file, header=0, sep="\t", dtype=dtype)

	i = 0
	while (i < n):
		r_perturb = perturb_df_translate(r_collection, i * step)
		yield t_collection, r_perturb
		i += 1


def perturb_reconsturction_scale(t_file, r_file, step=0.1, n=10):
	"""
	Perturb reconstruction file by scaling
	"""
	dtype = {'#chr': 'str',
			 'start': 'int',
			 'end': 'int',
			 'circ_id': 'str',
			 'estimated_cn': 'double',
			 'strand': 'str'}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)
	r_collection = pd.read_csv(r_file, header=0, sep="\t", dtype=dtype)

	i = 0
	while (i < n):
		r_perturb = perturb_df_scale(r_collection, i * step)
		yield t_collection, r_perturb
		i += 1


def perturb_reconsturction_extend_y(t_file, r_file, step=10, n=10):
	"""
	Perturb reconstruction file by scaling
	"""
	dtype = {'#chr': 'str',
			 'start': 'int',
			 'end': 'int',
			 'circ_id': 'str',
			 'estimated_cn': 'double',
			 'strand': 'str'}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)
	r_collection = pd.read_csv(r_file, header=0, sep="\t", dtype=dtype)

	i = 0
	while (i < n):
		r_perturb = perturb_df_extend_y(r_collection, i * step)
		yield t_collection, r_perturb
		i += 1


def perturb_reconsturction_rotate(t_file, r_file, step=10, n=10):
	"""
	Perturb reconstruction file by scaling
	"""
	dtype = {'#chr': 'str',
			 'start': 'int',
			 'end': 'int',
			 'circ_id': 'str',
			 'estimated_cn': 'double',
			 'strand': 'str'}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)
	r_collection = pd.read_csv(r_file, header=0, sep="\t", dtype=dtype)

	i = 0
	while (i < n):
		r_perturb = perturb_df_rotate(r_collection, i * step)
		yield t_collection, r_perturb
		i += 1


def perturb_reconsturction_rotate_2(t_file, r_file, step=10, grade=10, n=10):
	"""
	Perturb reconstruction file by scaling
	"""
	dtype = {'#chr': 'str',
			 'start': 'int',
			 'end': 'int',
			 'circ_id': 'str',
			 'estimated_cn': 'double',
			 'strand': 'str'}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)
	r_collection = pd.read_csv(r_file, header=0, sep="\t", dtype=dtype)

	i = 0
	while (i < n):
		r_perturb = perturb_df_rotate_2(r_collection, step=i * step, grade=i * grade)
		yield t_collection, r_perturb
		i += 1
