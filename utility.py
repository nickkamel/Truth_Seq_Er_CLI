from multiprocessing import Pool
import numpy as np

def multicore(tasks, func, num_processes, common_data):
    numTasksPerProcessFloor = len(tasks) // num_processes #This is rounded to nearest integer
    num_tasks_total_floor      = numTasksPerProcessFloor*num_processes
    tasksWhole              = tasks[0:num_tasks_total_floor]
    tasks_remainder          = tasks[num_tasks_total_floor:len(tasks)]

    results = list()
    if len(tasksWhole) !=0:
        n = numTasksPerProcessFloor
        processesTasks = [(tasksWhole[i:i + n], common_data) for i in range(0, num_tasks_total_floor, n)]
        pool = Pool(num_processes) 
        processes_results = pool.map(func, processesTasks)
        for p in processes_results:
            results += p
        pool.close()
        pool.join()

    if len(tasks_remainder) != 0:
        results += func([tasks_remainder, common_data])

    return results

def prefix_sum(input):
    accum = 0
    out = list()
    for i in input:
        accum += i
        out.append(accum)
    return out

def avg(inputList):
    return sum(inputList)/float(len(inputList))

#Sorts in ascending order 
def arg_sort(seq):
    return [i for (v, i) in sorted((v, i) for (i, v) in enumerate(seq))]

def flatten(l):
    flat_list = list()
    for sublist in l:
        for item in sublist:
            flat_list.append(item)
    return flat_list

def search_list_objs(targetList, attribute, targetValue):
    results = list()
    for storedObject in targetList:
        if getattr(storedObject, attribute) == targetValue:
            results.append(storedObject)
    return results 

def inv_list(targetList, value):
    results = [i for i, x in enumerate(targetList) if x == value]
    if results != []:
        return results[0]
    else:
        return None

#Takes a list of lists, where the lower level list represents the values for the multiple objectives
def non_dominated_sort(dataset):
    #print dataset
    numAttributes = len(dataset[0])
    list_size = len(dataset)
    num_ranked_individuals = 0
    dominated_counts = list()
    domination_sets = list()
    fronts = list()
    front0 = list()
    for i in range(list_size):
        domination_set = list()
        dominated_count = 0
        for j in range(list_size):
            numBetterAttributes_i = 0
            numBetterAttributes_j = 0
            numEqualAttributes = 0
            for attributeId in range(numAttributes):
                if dataset[i][attributeId] > dataset[j][attributeId]:
                    numBetterAttributes_i += 1
                elif dataset[i][attributeId] < dataset[j][attributeId]:
                    numBetterAttributes_j += 1
                else:
                    numEqualAttributes += 1
            #print dataset[i], dataset[j], numBetterAttributes_i, numBetterAttributes_j, numEqualAttributes
            
            if numBetterAttributes_j == 0 and numEqualAttributes != numAttributes: #if j is never better than i and i is better than j at at least 1 fitness
                 domination_set.append(j)
            elif numBetterAttributes_i == 0 and numEqualAttributes != numAttributes: #j dominates i
                 dominated_count += 1
                        
        if dominated_count == 0:
            front0.append(i)
            num_ranked_individuals += 1
        domination_sets.append(domination_set)
        dominated_counts.append(dominated_count)

    fronts.append(front0)
    front_counter = 0
    while(num_ranked_individuals != list_size):
        new_front = list()
        for i in fronts[front_counter]:
            for j in domination_sets[i]:
                dominated_counts[j] -= 1
                if dominated_counts[j] == 0:
                    new_front.append(j)
                    num_ranked_individuals +=1                        
        fronts.append(new_front)
        front_counter +=1

    return fronts

def dot_distance(v1, v2):
    return np.sum(np.abs(np.subtract(v1, v2)) )

def dot_distance_batch(batch):
    tasks, common_data = batch
    results = list()
    for task in tasks:
        results.append(dot_distance(*task))
    return results

def remove_from_list(input_list, items):
    for item in items:
        if item in input_list:
            input_list.remove(item)
    return input_list

def bin_to_dec(binary):
    return sum(val*(2**idx) for idx, val in enumerate(reversed(binary)))