from scipy.io import loadmat
graphData = loadmat('MulNet.mat')
#the node are 1 to 6407 and not 0 to 6406 and they are: yt-tw-ff

#proof that friendfeed is directed: https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.researchgate.net%2Ffigure%2FComparative-performances-on-FriendFeed-damp-ing-factor-depending-on-the-average-path_tbl9_265729958&psig=AOvVaw1QTIc9sX5GrBFkzqzobZPl&ust=1646661336214000&source=images&cd=vfe&ved=0CAsQjRxqFwoTCMjBntfRsfYCFQAAAAAdAAAAABAD

print(graphData.keys())
print(graphData['NodeDegree'])
'''
for i in range(0, 6407):
	for j in range(0, 3):
	print(graphData['ListGraph'][i][j])

print(graphData['NodesFeach'])'''

# now create three files tw.txt, ff.txt and yt.txt:
#yt:
'''
for i in range(0, 6407):
	currentNodeFriendsList = list(graphData['ListGraph'][i][0])[0]

	for friend in range(0, len(currentNodeFriendsList)):
		with open('yt.txt', 'a') as file:
			file.write(str(i+1) + " " + str(currentNodeFriendsList[friend]) + "\n")

#tw:

for i in range(0, 6407):
	currentNodeFriendsList = list(graphData['ListGraph'][i][1])[0]

	for friend in range(0, len(currentNodeFriendsList)):
		with open('tw.txt', 'a') as file:
			file.write(str(i+1) + " " + str(currentNodeFriendsList[friend]) + "\n")

#ff:

for i in range(0, 6407):
	currentNodeFriendsList = list(graphData['ListGraph'][i][2])[0]

	for friend in range(0, len(currentNodeFriendsList)):
		with open('ff.txt', 'a') as file:
			file.write(str(i+1) + " " + str(currentNodeFriendsList[friend]) + "\n")
'''
