
import numpy as np
import pandas as pd
 
def min_dist_calc(x1,y1,z1,x2,y2,z2,box_size):
    dist1=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    #|1|2|3|
    #|4|x|5|
    #|6|7|8|
    #box 4,5
    if x1>x2:
        dist2=np.sqrt((x1-x2-box_size)**2+(y1-y2)**2+(z1-z2)**2)
    else:
        dist2=np.sqrt((x1-x2+box_size)**2+(y1-y2)**2+(z1-z2)**2)
    #|1|2|3|
    #|4|x|5|
    #|6|7|8|
    #box 2,7
    if y1>y2:
        dist3=np.sqrt((x1-x2)**2+(y1-y2-box_size)**2+(z1-z2)**2)
    else:
        dist3=np.sqrt((x1-x2)**2+(y1-y2+box_size)**2+(z1-z2)**2)
    #|1|2|3|
    #|4|x|5|
    #|6|7|8|
    #box 3,6
    if x1>x2 and y1>y2:
        dist4=np.sqrt((x1-x2-box_size)**2+(y1-y2-box_size)**2+(z1-z2)**2)
    else:
        dist4=np.sqrt((x1-x2+box_size)**2+(y1-y2+box_size)**2+(z1-z2)**2)
    #|1|2|3|
    #|4|x|5|
    #|6|7|8|
    #box 1,8
    if x1<x2 and y1>y2:
        dist5=np.sqrt((x1-x2+box_size)**2+(y1-y2-box_size)**2+(z1-z2)**2)
    else:
        dist5=np.sqrt((x1-x2-box_size)**2+(y1-y2+box_size)**2+(z1-z2)**2)
    min_dist=min(dist1,dist2,dist3,dist4,dist5)
    return min_dist

def dist_form_3d(x1,y1,z1,x2,y2,z2):
    dist=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    return dist

def dist_form_2d(y1,z1,y2,z2):
    dist=np.sqrt((y1-y2)**2+(z1-z2)**2)
    return dist

def center_calc(array):
    try:
        if len(array)%3!=0:
            print("EXCEPTION:::missing coordinate")
        else:
            num_coord=len(array)/3
    except:
        print("EXCEPTION:::format should be array [x1,y1,z1, ... xn,yn,zn]")
    x=0
    y=0
    z=0
    for cnt in np.arange(num_coord):
        xcnt=cnt*3
        ycnt=cnt*3+1
        zcnt=cnt*3+2
        x=x+array[int(xcnt)]
        y=y+array[int(ycnt)]
        z=z+array[int(zcnt)]
    x_cen=x/num_coord
    y_cen=y/num_coord
    z_cen=z/num_coord
    center=[x_cen,y_cen,z_cen]
    return center

def calc_rel_arm_angle(b,angleNum,armNum,protNum):
    ang=pd.DataFrame()
    xCoord=[]
    yCoord=[]
    ##################### ANGLE CALCULATION PREPARATION
    ############### X coord zeroed to origin
    for protCnt in np.arange(protNum):
        protXColumn=protCnt*3+1
        for armCnt in np.arange(armNum):
            for angleCnt in np.arange(angleNum):
                angleXColumn=protNum*3+protCnt*armNum*angleNum*3+armCnt*angleNum*3+angleCnt*3+1
                xCoord.append(b[angleXColumn]-b[protXColumn])
    ############### Y coord zeroed to origin
    for protCnt in np.arange(protNum):
        protYColumn=protCnt*3+2
        for armCnt in np.arange(armNum):
            for angleCnt in np.arange(angleNum):
                angleYColumn=protNum*3+protCnt*armNum*angleNum*3+armCnt*angleNum*3+angleCnt*3+2
                yCoord.append(b[angleYColumn]-b[protYColumn])
    ##################### ANGLE CALCULATION
    frameAngles=np.arctan2(yCoord,xCoord)*180/np.pi
    ################## NUMBER OF COMPARISONS DEPENDENT ON NUMBER OF PROTEINS
    comparisonNum=0 #DEFAULT SETTING DO NOT CHANGE
    comparison_not_calculated=True
    if comparison_not_calculated:
        for protCnt in np.arange(protNum):
            comparisonNum += protCnt
        comparison_not_calculated=False
    ################## FINDING THE MINIMUM ANGLE WITH DIFFERENT ANGLE CALCULATIONS AND CALCULATE THEIR AVERAGE
    relative_angle=[]
    angle_type=[]
    for angleCnt in np.arange(angleNum):
        prot_prot_angle_mean=[]
        for i,protCnt in enumerate(np.arange(protNum)):
            for f,otherProtCnt in enumerate(np.arange(protNum)):
                if (i != f) and (i > f):
                    compare_min_list=[]
                    for armCnt in np.arange(armNum):
                        comparison_list=[]
                        for otherProtArmCnt in np.arange(armNum):
                            angle1=frameAngles[angleCnt+(armCnt*angleNum)+(protCnt*angleNum*armNum)]
                            angle2=frameAngles[(otherProtCnt*angleNum*armNum)+(otherProtArmCnt)*angleNum+angleCnt]
                            comparison_list.append(abs(angle1-angle2))
                            comparison_list.append(360-abs(angle1-angle2))
                        compare_min_list.append(min(comparison_list)) # compare 1 to other protein arms and pick the relative minimum
                        del comparison_list
                    prot_prot_angle_mean.append(np.mean(compare_min_list)) # compare the relative angle calculation of both proteins using all arms
                    del compare_min_list
        angle_type.append(prot_prot_angle_mean) # saving the calculation for 1 angle type (tip or middle of arm)
        del prot_prot_angle_mean
    for comparison in np.arange(comparisonNum):
        for angleCnt in np.arange(angleNum-1):
            relative_angle.append((angle_type[angleCnt][comparison]+angle_type[angleCnt+1][comparison])/angleNum)
    #####################
    ang = ang.append(
        {
            'Frame':b[0],
            'rel_ang':relative_angle[0],
        },
        ignore_index=True
    )
    return ang

def rolling_average(data, window_size):
    # Create a 1D window of ones to compute the average
    window = np.ones(window_size) / window_size
    # Use the 'valid' mode to ignore boundary effects
    return np.convolve(data, window, mode='valid')

def rolling_average_by_simname(df, window_size, ignore_index=False):
    # Group the DataFrame by 'simName' column
    grouped = df.groupby('simName')
    
    # Initialize a list to store rolling average DataFrames
    rolling_avg_dfs = []
    
    # Iterate through each group and perform rolling average
    for name, group in grouped:
        data = group.drop('simName', axis=1)  # Remove 'simName' column
        rolling_avg_data = data.rolling(window=window_size, axis=0).mean()
        
        # Add 'simName' back to the rolling average data
        rolling_avg_data['simName'] = name
        
        # Append to the list of rolling_avg_dfs
        rolling_avg_dfs.append(rolling_avg_data)
    
    # Concatenate all rolling average DataFrames into one
    rolling_avg_df = pd.concat(rolling_avg_dfs, ignore_index=ignore_index)
    
    return rolling_avg_df

#def rolling_average_by_simname(df, window_size, ignore_index=False):
#    # Define the rolling average function
#    def rolling_average(col, window_size):
#        return col.rolling(window=window_size, min_periods=1).mean()
#
#    # Group the DataFrame by 'simName' column
#    grouped = df.groupby('simName')
#
#    # Initialize a new DataFrame to store rolling averages
#    rolling_avg_df = pd.DataFrame(columns=df.columns)
#
#    # Iterate through each group and perform rolling average
#    for name, group in grouped:
#        data = group.drop('simName', axis=1)  # Remove 'simName' column
#        rolling_avg_data = data.apply(lambda col: rolling_average(col, window_size), axis=0)
#
#        # Add 'simName' back to the rolling average data
#        rolling_avg_data['simName'] = name
#
#        # Append to the rolling_avg_df
#        rolling_avg_df = rolling_avg_df.append(rolling_avg_data, ignore_index=ignore_index)
#
#    return rolling_avg_df
