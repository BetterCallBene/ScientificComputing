import csv
import sys
import numpy as np
import matplotlib.pyplot as plt

def csvFile(filename):
	data = []
	with open(filename, 'rb') as csvFile:
		reader = csv.reader(csvFile, delimiter=',')
		for row in reader:
			l = []
			for elm in row:
				l.append(float(elm))
			data.append(l)

	arr = np.array(data, float)
	return arr

def plotData(data):
	arr = np.transpose(data)

	plt.figure(1)
	plt.subplot(311)
	plt.plot(arr[0], arr[1], 'ro')
	plt.subplot(312)
	plt.plot(arr[0], arr[2], 'b--')
	plt.subplot(313)
	plt.plot(arr[0], arr[3], 'gs')
	plt.show()



if __name__ == "__main__":
	if sys.argv[1] is not None:
		filename = sys.argv[1]
		data = csvFile(filename)
		plotData(data)
	else:
		print 'Bitte Dateiname angeben'