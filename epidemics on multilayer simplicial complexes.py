#!C:\pythonCode
# -*- coding: utf-8 -*-
# @Author : fjf
# @Software: PyCharm
# @Profile :
# ******constant variables******
# ***********************
ITS = 0.4  # μ:disease recovery rate
ATU = 0.8  # δ:information loss rate
lambda_delta = 0.6  # λ_δ:control the information spreading rate of 2-simplex
# ***********************
experiment_times = 1
PC = 1
# ***********************
experiment_date = "2022021601"
file_network_triangles = "your route"
file_network_down = 'your route'
file_network_up = "your route"
file_path = "your route"
result_path = "your route" + experiment_date + "_"
# ***********************
N = 1000
I_initial = 0.01
T_step = 100000
A1_step = 50
A1_init = round(2/100, 2)
A2_step = 50
A2_init = round(2/100, 2)
REDUCE_A1 = 0  # work as γ,β(A)= γβ(U)
T = 20

# ******open file******
fA = open(result_path + "sis_mmca_evolution_a_" + str(experiment_times) + ".txt", "w")
fI = open(result_path + "sis_mmca_evolution_i_" + str(experiment_times) + ".txt", "w")

# ******write the data of networks to the 2-dimension list******
f_network_up = open(file_path + file_network_up)
f_network_triangles = open(file_path + file_network_triangles)
f_network_down = open(file_path + file_network_down)
list_network_up = []
triangles_list = []
list_network_down = []
for line in f_network_up:
    temp = line.split()
    temp = [int(x) for x in temp]
    list_network_up.append(temp)
for line in f_network_triangles:
    temp = line.split()
    temp = [int(x) for x in temp]
    triangles_list.append(temp)
for line in f_network_down:
    temp = line.split()
    temp = [int(x) for x in temp]
    list_network_down.append(temp)
print("2-dimension list is finished")

# ******construct the adjacency matrix******
matrix_network_up = []
matrix_network_down = []
for i in range(N):
    matrix_network_up.append([0 for n in range(N)])
    matrix_network_down.append([0 for n in range(N)])
for i in range(len(list_network_up)):
    m = list_network_up[i][0]
    n = list_network_up[i][1]
    matrix_network_up[m][n] = 1
    matrix_network_up[n][m] = 1
for i in range(len(list_network_down)):
    m = list_network_down[i][0]
    n = list_network_down[i][1]
    matrix_network_down[m][n] = 1
    matrix_network_down[n][m] = 1
print("the adjacency matrix is finished")

# ******construct the adjacency list******
adjacency_list_up = []
adjacency_list_down = []
for i in range(N):
    temp_up = []
    temp_down = []
    for j in range(N):
        if matrix_network_up[i][j] == 1:
            temp_up.append(j)
        else:
            pass
        if matrix_network_down[i][j] == 1:
            temp_down.append(j)
        else:
            pass
    adjacency_list_up.append(temp_up)
    adjacency_list_down.append(temp_down)
print("the adjacency list is finished")

# ******construct the adjacency list of 2-simplex******
adjacency_triangles_list = [[] for i in range(N)]
for i in range(len(triangles_list)):
    sublist = triangles_list[i]
    for j in range(3):
        if j == 0:
            cut_list = list()
            cut_list.append(sublist[1])
            cut_list.append(sublist[2])
            adjacency_triangles_list[sublist[0]].append(cut_list)
        elif j == 1:
            cut_list = list()
            cut_list.append(sublist[0])
            cut_list.append(sublist[2])
            adjacency_triangles_list[sublist[1]].append(cut_list)
        elif j == 2:
            cut_list = list()
            cut_list.append(sublist[0])
            cut_list.append(sublist[1])
            adjacency_triangles_list[sublist[2]].append(cut_list)
        else:
            pass

# ******start MMCA circle******
for i in range(A1_step+1):
    UTA = round(i*A1_init, 4)  # λ:information spreading rate
    for j in range(A2_step+1):
        STA_InfectionStrength = round(A2_init * j, 4)  # β:disease transmission rate
        A_InfectionStrength = round(STA_InfectionStrength * REDUCE_A1, 4)
        # init the proportion of all states
        US_Size = [(1 - I_initial) for n in range(N)]
        AS_Size = [0 for n in range(N)]
        AI_Size = [I_initial for n in range(N)]
        A_Size = [I_initial for n in range(N)]
        A_nodes_percent_t = 0
        I_nodes_percent_t = 0
        stop = 0
        print("UTA：" + str(UTA) + "，STI：" + str(STA_InfectionStrength))

        for t in range(T_step):
            A_Count = 0
            I_Count = 0
            A_size_pre = A_Size.copy()
            AI_size_pre = AI_Size.copy()
            for n in range(N):
                # r_i^1 (t):the possibility of nodes not be informed by aware neighbours
                r_A_neighbour_temp = 1
                for neighbor in adjacency_list_up[n]:
                    r_A_neighbour_temp = (1 - A_size_pre[neighbor] * UTA) * r_A_neighbour_temp
                # r_i^2 (t):the possibility of nodes not be informed by 2-simplex neighbours
                r_A_triangles_temp = 1
                real_kD_up = 3. * len(triangles_list) / len(adjacency_list_up)
                # λ^*: the information spreading rate of 2-simplex
                lambda_u = lambda_delta * ATU / real_kD_up
                for s_list in adjacency_triangles_list[n]:
                    r_A_triangles_temp = (1 - A_size_pre[s_list[0]] * A_size_pre[
                        s_list[1]] * lambda_u) * r_A_triangles_temp
                # r_i (t):the possibility of nodes not be informed
                u_per = r_A_neighbour_temp * r_A_triangles_temp
                A_per = 1 - u_per
                # q_i^A (t):the possibility of aware nodes not be infected
                # q_i^U (t):the possibility of unaware nodes not be infected
                q_A = 1
                q_U = 1
                for neighbor in adjacency_list_down[n]:
                    q_A = (1 - AI_size_pre[neighbor] * A_InfectionStrength) * q_A
                    q_U = (1 - AI_size_pre[neighbor] * STA_InfectionStrength) * q_U

                temp_AI, temp_US, temp_AS = AI_Size[n], US_Size[n], AS_Size[n]
                US_Size[n] = temp_US * u_per * q_U + temp_AI * ATU * ITS + temp_AS * ATU * q_U
                AI_Size[n] = temp_US * A_per * (1 - q_A) + temp_US * u_per * (1 - q_U) + temp_AI * (
                        1 - ITS) + temp_AS * ATU * (1 - q_U) + temp_AS * (1 - ATU) * (1 - q_A)
                AS_Size[n] = temp_US * A_per * q_A + temp_AI * (1 - ATU) * ITS + temp_AS * (1 - ATU) * q_A
                A_Size[n] = AS_Size[n] + AI_Size[n]
                # %%%%%%check%%%%%%
                checksum = AI_Size[n] + US_Size[n] + AS_Size[n]
                if checksum > 1.00001 or checksum < 0.99999:
                    print(AI_Size[n], US_Size[n], AS_Size[n])
                    print(AI_Size[n] + US_Size[n] + AS_Size[n])
                    input("check errors！")
                else:
                    pass

                A_Count += A_Size[n]
                I_Count += AI_Size[n]
            A_size_step = A_Count / N
            I_size_step = I_Count / N
            # ******judge if the numbers of infectious nodes are steady or not******
            if t >= 300 and stop != 1:
                if t % 30 == 0:
                    stop_list = []
                    I_size_first = I_size_step
                else:
                    pass
                stop_list.append(I_size_step)
                if t % 30 == 14:
                    I_size_second = I_size_step
                else:
                    pass
                if t % 30 == 29:
                    diff_sum = 0
                    I_size_third = I_size_step
                    for st in range(len(stop_list) - 1):
                        diff = stop_list[st] - stop_list[st + 1]
                        diff_sum += diff
                    if abs(I_size_first - I_size_third) < 0.005 and abs(I_size_first - I_size_second) < 0.005 and abs(
                            I_size_second - I_size_third) < 0.005 and abs(diff_sum) < 0.005:
                        stop_t = t
                        stop = 1
                    else:
                        pass
                else:
                    pass
            else:
                pass
            # ******if steady,  stop circles******
            if stop == 1:
                if t < stop_t + T:
                    A_nodes_percent_t += A_size_step
                    I_nodes_percent_t += I_size_step
                else:
                    break
            else:
                pass
            if t >= T_step - T:
                input("beyond the maximum time scale！")
            else:
                pass

        # ******write the result to file******
        str_A = str(UTA) + '  ' + str(STA_InfectionStrength) + "    " + str(round(A_nodes_percent_t / T, 4)) + "\n"
        str_I = str(UTA) + '  ' + str(STA_InfectionStrength) + "    " + str(round(I_nodes_percent_t / T, 4)) + "\n"
        fA.write(str_A)
        fI.write(str_I)

# ******close file******
fA.close()
fI.close()
f_network_up.close()
f_network_down.close()
f_network_triangles.close()


