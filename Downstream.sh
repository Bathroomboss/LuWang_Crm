import sys
# 输入TE 和 dm6 的 paf文件， 然后输出4个文件, partial reads，spanning reads，spanning reads 的组成成分
# 比之前的版本多了一个去除‘translocation’的功能

# 修正了以前spanning reads寻找错误的bug
# 删掉一条reads多个片段（>3）回贴到相同的TE的信息,完善了translocation 对 genome位置的剔除
# 并且保留 overlap中较长的片段
# 考虑完全包含，如果某个片段完全包含另一个片段，则较长的是真实的，如果不是完全包含，overlap过多的将被删除
# 在print partial insertion 的时候，也需要print出 TE insertion 的位点

TE = open(sys.argv[1],"r")
TE_List = []
for line in TE.readlines():
    line = line.replace('\t',',').replace("\n","")
    #print(line)
    TElinelist = line.split(',')
    new_TElinelist = [int(n) if n.isdigit() else n for n in TElinelist]
    TE_List.append(new_TElinelist)

dm6 = open(sys.argv[2],"r")
dm6_List = []
for line in dm6.readlines():
    line = line.replace('\t',',').replace("\n","")
    #print(line)
    dm6linelist = line.split(',')
    new_dm6linelist = [int(n) if n.isdigit() else n for n in dm6linelist]
    dm6_List.append(new_dm6linelist)
# print(dm6_List)
# 将TE的list整成字典
te_dict = {}
for te_sublist in TE_List:
    te_read_name = te_sublist[0]
    if te_read_name in te_dict:
        te_dict[te_read_name].append(te_sublist)
    else:
        te_dict[te_read_name] = [te_sublist]

# 删掉一条reads多个片段比对到同一个TE的情况
new_te_dict = {}
for te_read_name, te_sublists in te_dict.items():
    te_name_counts = {}
    for sublist in te_sublists:
        te_name = sublist[5]
        if te_name in te_name_counts:
            te_name_counts[te_name] += 1
        else:
            te_name_counts[te_name] = 1
    filtered_value = [sublist for sublist in te_sublists if te_name_counts[sublist[5]] <= 3]
    new_te_dict[te_read_name] = filtered_value
te_dict = new_te_dict
# 得到过滤之后的te paf 文件
new_TE_paf = list(te_dict.values())

# # 检测‘translocation’
# 如果一个reads 的一段map到了TE，但是这一段的一部分也是基因组的某个位置，则不认为是insertion，将其删除
# ”一部分定义为“，相互 overlap的区间，不得超过 reads map到TE长度的50%
def compute_overlap_length(a, b, x, y):
    # a,b 为 TE 在reads 的区间
    # x,y 为 dm6 在reads 的区间
    if b < x or y < a:
        return 0  # 区间a, b和区间x, y没有重叠
    else:
        start = max(a, x)
        end = min(b, y)
        overlap_length = end - start
        return overlap_length

te_sublists_to_remove = []
dm6_sublists_to_remove = []
#te_overlap_info = []
for dm6_sublist in dm6_List:
    dm6_read_name = dm6_sublist[0]
    if dm6_read_name in te_dict:
        te_sublists = te_dict[dm6_read_name]
        #te_sublists_to_remove.extend([te_sublist for te_sublist in te_sublists if compute_overlap_length(te_sublist[2], te_sublist[3], dm6_sublist[2], dm6_sublist[3]) >= 0.5 * (te_sublist[3] - te_sublist[2])])
        for te_sublist in te_sublists:
            overlap_L = compute_overlap_length(te_sublist[2], te_sublist[3], dm6_sublist[2], dm6_sublist[3]) 
            if overlap_L != 0:
                dm6_L = dm6_sublist[3] - dm6_sublist[2]
                TE_L = te_sublist[3] - te_sublist[2]
                dm6_ratio = overlap_L / dm6_L
                TE_ratio = overlap_L / TE_L
                if TE_ratio == 1 and dm6_ratio == 1: # 删掉完全一致的片段
                    te_sublists_to_remove.append(te_sublist) 
                    dm6_sublists_to_remove.append(dm6_sublist)
                if TE_ratio >= 0.5 and not (te_sublist[2] <= dm6_sublist[2] and te_sublist[3] >= dm6_sublist[3]): # 如果TE片段 overlap过多，并且不是完全包含dm6片段
                    te_sublists_to_remove.append(te_sublist)
                if dm6_ratio >= 0.5 and not (dm6_sublist[2] <= te_sublist[2] and dm6_sublist[3] >= te_sublist[3]): # 如果dm6片段 overlap过多，并且不是完全包含TE片段
                    dm6_sublists_to_remove.append(dm6_sublist)

# 从 te_sublists 中删除需要删除的 te_sublist
for dm6_read_name in te_dict:
    te_sublists = te_dict[dm6_read_name]
    #for te_sublist in te_sublists_to_remove:
    te_sublists = [te_sublist for te_sublist in te_sublists if te_sublist not in te_sublists_to_remove]
    te_dict[dm6_read_name] = te_sublists

# 删除不需要的dm6_list
dm6_List = [sublist for sublist in dm6_List if sublist not in dm6_sublists_to_remove]

# 根据reads的名称进行比较
TE_left_list = []
TE_right_list = []
for dm6_sublist in dm6_List:
    dm6_read_name = dm6_sublist[0]
    if dm6_read_name in te_dict:
        te_sublists = te_dict[dm6_read_name]
        for te_sublist in te_sublists:
            if -20 < te_sublist[2] - dm6_sublist[3] < 100: # TE在右边
                #print(dm6_sublist[0],dm6_sublist[2],dm6_sublist[3],te_sublist[2],te_sublist[3],
                #      dm6_sublist[5],dm6_sublist[7],dm6_sublist[8],dm6_sublist[4], te_sublist[5],te_sublist[7],te_sublist[8],te_sublist[4])
                TE_right_list.append([dm6_sublist[0],dm6_sublist[2],dm6_sublist[3],te_sublist[2],te_sublist[3],
                      dm6_sublist[5],dm6_sublist[7],dm6_sublist[8],dm6_sublist[4], te_sublist[5],te_sublist[7],te_sublist[8],te_sublist[4]])
            elif -20 < dm6_sublist[2] - te_sublist[3] < 100: # TE在左边
                #print(dm6_sublist[0], te_sublist[2], te_sublist[3], dm6_sublist[2], dm6_sublist[3],
                #      te_sublist[5], te_sublist[7], te_sublist[8], te_sublist[4],dm6_sublist[5], dm6_sublist[7], dm6_sublist[8], dm6_sublist[4])
                TE_left_list.append([dm6_sublist[0], te_sublist[2], te_sublist[3], dm6_sublist[2], dm6_sublist[3],
                      te_sublist[5], te_sublist[7], te_sublist[8], te_sublist[4],dm6_sublist[5], dm6_sublist[7], dm6_sublist[8], dm6_sublist[4]])
#print(TE_left_list)
#print(TE_right_list)
# 对于left reads 和 right reads 查找 spanning Insertion
TE_left_dict = {}
for te_left_sublist in TE_left_list:
    te_left_reads = te_left_sublist[0]
    if te_left_reads in TE_left_dict:
        TE_left_dict[te_left_reads].append(te_left_sublist)
    else:
        TE_left_dict[te_left_reads] = [te_left_sublist]
#print(TE_left_dict)
spanning_insertion = []
Spanning_Constituent_element = []
for te_right_sublist in TE_right_list:
    te_right_reads = te_right_sublist[0]
    if te_right_reads in TE_left_dict:
        te_left_sublists = TE_left_dict[te_right_reads]
        for te_left_sublist in te_left_sublists:
            if te_left_sublist[5:9] == te_right_sublist[9:13] and te_right_sublist[3:5] == te_left_sublist[1:3] and te_left_sublist[9] == te_right_sublist[5] and te_left_sublist[12] == te_right_sublist[8]:
                if te_left_sublist[12] == "+":
                    if abs(te_left_sublist[10] - te_right_sublist[7]) < 200:
                        spanning_insertion.append([te_right_sublist[0],te_left_sublist[9],min(te_left_sublist[10],te_right_sublist[7]),max(te_left_sublist[10],te_right_sublist[7]),te_left_sublist[12],te_left_sublist[5],te_left_sublist[6],te_left_sublist[7],te_left_sublist[8]])
                        Spanning_Constituent_element.append((te_right_sublist, te_left_sublist))
                elif te_left_sublist[12] == "-":
                    if abs(te_right_sublist[6] - te_left_sublist[11]) < 200:
                        spanning_insertion.append([te_right_sublist[0],te_left_sublist[9],min(te_right_sublist[6] ,te_left_sublist[11]),max(te_right_sublist[6],te_left_sublist[11]),te_left_sublist[12],te_left_sublist[5],te_left_sublist[6],te_left_sublist[7],te_left_sublist[8]])
                        Spanning_Constituent_element.append((te_right_sublist, te_left_sublist))

#print(Spanning_Constituent_element)
for right_sublist, left_sublist in Spanning_Constituent_element:
    TE_right_list.remove(right_sublist)
    TE_left_list.remove(left_sublist)

with open(sys.argv[3], "w") as file:
    file.write("TE_right_list:\n")
    for sublist in TE_right_list:
        file.write(str(sublist) + "\n")
    file.write("\nTE_left_list:\n")
    for sublist in TE_left_list:
        file.write(str(sublist) + "\n")
# 存储 spanning_insertion
with open(sys.argv[4], "w") as file:
    for sublist in spanning_insertion:
        file.write(str(sublist) + "\n")
# 存储 Spanning_Constituent_element
with open(sys.argv[5], "w") as file:
    for right_sublist, left_sublist in Spanning_Constituent_element:
        file.write("Right Sublist: " + str(right_sublist) + "\n")
        file.write("Left Sublist: " + str(left_sublist) + "\n")
        file.write("\n")

#with open(sys.argv[6], "w") as file:
#    for reads in keys_to_delete:
#        file.write(str(reads) + "\n")

#with open(sys.argv[6], "w") as file:
#    for reads, value in zip(keys_to_delete, TE_circle_detect):
#        file.write(str(reads) + " " + str(value) + "\n")


Sum_partial_list = []
for sublist in TE_right_list:
    if sublist[9] == "HMS-Beagle-rpl14-GFPreporter:Gypsy":
        if sublist[10] <= 8546 and sublist[11] >= 6201:
            sublist[9] = "HMS-Beagle-rpl14-GFPreporter:Gypsy"
        else:
            sublist[9] = "HMS-Beagle"
    if sublist[8] == "+":
        dm6_insert_site = sublist[7]
        if sublist[12] == "+":
            TE_insert_site = sublist[10]
        elif sublist[12] == "-":
            TE_insert_site = sublist[11]
        extracted_elements = [sublist[5],(dm6_insert_site//200)*200,(dm6_insert_site//200)*200+200,sublist[8], sublist[9],dm6_insert_site,TE_insert_site,sublist[0]]
        Sum_partial_list.append(extracted_elements)
    elif sublist[8] == "-":
        dm6_insert_site = sublist[6]
        if sublist[12] == "+":
            TE_insert_site = sublist[10]
        elif sublist[12] == "-":
            TE_insert_site = sublist[11]
        extracted_elements = [sublist[5],(dm6_insert_site//200)*200,(dm6_insert_site//200)*200+200, sublist[8], sublist[9],dm6_insert_site,TE_insert_site,sublist[0]]
        Sum_partial_list.append(extracted_elements)
# 对于 TE_left_list 中的每个小列表
for sublist in TE_left_list:
    if sublist[5] == "HMS-Beagle-rpl14-GFPreporter:Gypsy":
        if sublist[6] <= 8546 and sublist[7] >= 6201:
            sublist[5] = "HMS-Beagle-rpl14-GFPreporter:Gypsy"
        else:
            sublist[5] = "HMS-Beagle"
    if sublist[12] == "+":
        dm6_insert_site = sublist[10]
        if sublist[8] == "+":
            TE_insert_site = sublist[7]
        elif sublist[8] == "-":
            TE_insert_site = sublist[6]
        extracted_elements = [sublist[9], (dm6_insert_site//200)*200,(dm6_insert_site//200)*200+200, sublist[12], sublist[5],dm6_insert_site,TE_insert_site,sublist[0]]
        Sum_partial_list.append(extracted_elements)
    elif sublist[12] == "-":
        dm6_insert_site = sublist[11]
        if sublist[8] == "+":
            TE_insert_site = sublist[7]
        elif sublist[8] == "-":
            TE_insert_site = sublist[6]
        extracted_elements = [sublist[9], (dm6_insert_site//200)*200,(dm6_insert_site//200)*200+200, sublist[12], sublist[5],dm6_insert_site,TE_insert_site,sublist[0]]
        Sum_partial_list.append(extracted_elements)

with open(sys.argv[6], "w") as file:
    for sums in Sum_partial_list:
        file.write(str(sums) + "\n")

with open(sys.argv[7], "w") as file:
    for sublist in new_TE_paf:
        file.write(str(sublist) + "\n")

#with open(sys.argv[8], "w") as file:
#    for rm in  te_sublists_to_remove:
#        file.write(str(rm) + "\n")
#with open(sys.argv[9], "w") as file:
#    for overlap_info in te_overlap_info:
#        file.write(str(overlap_info) + "\n")

print("Done! Good Job!")
