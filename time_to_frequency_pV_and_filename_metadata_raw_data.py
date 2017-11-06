import os
import glob
import pandas
import numpy as np
from scipy import signal
from datetime import datetime

# raw csv or modified (difference in line number start of file, header_row)
raw_csv = True

# the date that MEA2 was seeded
date_seeded = "2016-10-11"

# simple low pass filter
#frequency_treshold = len(f) # [Hz]
frequency_treshold = 5000 # [Hz]

# get filenames
raw_filenames = sorted(glob.glob('*.csv'))

# sampling info
Fs = 10000 # [1/s] # 10 kHz
dt = 1/float(Fs) # [s]

# create output directory PSD if it doesn't exist
if not os.path.exists("PSD"):
    os.makedirs("PSD")

# define functions
def find_index_of_list_containing_string(l, s):
    for e in l:
        if s in e:
            return l.index(e)

def extract_experiment_stimulation(file_to_import):   
    start_string_type_1 = "MEA 2 "
    start_string_type_2 = "MEA2 "
    
    if start_string_type_1 in file_to_import:
        index_start = file_to_import.index(start_string_type_1) + 6
    elif start_string_type_2 in file_to_import:
        index_start = file_to_import.index(start_string_type_2) + 5
    else:
        raise ValueError('could not find start string patterns in %s' % (file_to_import))
        
    stop_string = "Recording-0"
    
    if not stop_string in file_to_import:
        raise ValueError('could not find stop string pattern in %s' % (file_to_import))
    
    index_stop = file_to_import.index(stop_string) - 1
    
    return file_to_import[index_start:index_stop]

def lies_in_frequency_band(Hz):
    # return conventional frequency band for Hz
    # according to neuro science:
    
    # delta: [0.5, 4>
    # theta: [4, 8>
    # alpha: [8, 13>
    # beta: [13, 30>
    # gamma: [30, 120>
    # unspecified: [0, 0.5> or [120, inf>
    
    # assuming appropriate preprocessing
    
    if Hz >= 0.5 and Hz < 4:
        return "delta"
    elif Hz >= 4 and Hz < 8:
        return "theta"
    elif Hz >= 8 and Hz < 13:
        return "alpha"
    elif Hz >= 13 and Hz < 30:
        return "beta"
    elif Hz >= 30 and Hz < 120:
        return "gamma"
    elif Hz >= 0 and (Hz >= 120 or Hz < 0.5):
        return "unspecified"
    else:
        raise ValueError("could not identify frequency band of %i" % Hz)

def extract_experiment_date(file_to_import):
    try:
        index_start = file_to_import.index("/") + 1
    except ValueError:
        index_start = 0
    return file_to_import[index_start:index_start+10]

def days_between(d1, d2):
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days) + 1 # regard one day experiment as one day long

def get_total_energy_for_frequency(P, Hz, Fs, time):
    # total energy for a specific frequency
    # assume Hz is an integer
    return P[Hz]*Fs*time

def make_unscrambler_friendly_list(numpy_array_or_list):
    # converts number into string and replaces . comma with , comma
    return [ str(number).replace('.',',') for number in numpy_array_or_list ]

def calculate_signal_duration(L, dt):
    # return signal duration in seconds using
    # L: the number of samples in the singal S. len(S)
    # dt: the periodic time of the sampling. 1/Fs
    return L*dt

def calculate_signal_mean(S):
    return np.mean(S)
    
def calculate_signal_min(S):
    return np.min(S)

def calculate_signal_max(S):
    return np.max(S)

def calculate_total_signal_energy(P_downsampled, f_downsampled, Fs, signal_duration):
    # total signal energy as a summation of the total
    # signal energy for each frequencies
    # using get_total_energy_for_frequency(P, Hz, Fs, signal_duration)

    # Fs: sampling rate [1/s]
    # signal duration: [s]
    return np.sum([get_total_energy_for_frequency(P_downsampled, int(np.ceil(Hz)), Fs, signal_duration) for Hz in f_downsampled])




# --- outer loop start ---

#for experiment_index in range(0, 1):#len(raw_filenames)):
for experiment_index in range(0, len(raw_filenames)):

    # select experiment index
    #experiment_index = 0

    # make empty pandas DataFrame to later save to csv file
    P_pandas_total_df = pandas.DataFrame()

    # chose file to import
    file_to_import = raw_filenames[experiment_index]
    print "Experiment info: Select experiment " + file_to_import

    if raw_csv:
        header_row = 5
    elif not raw_csv: # modified 1s samples A-E (TTK19)
        header_row = 0

    # Load the 61 column file
    print "Experiment info: Load the experiment" 
    raw_file = pandas.read_csv(file_to_import, sep=",", header=header_row, low_memory=True)

    # Extract header. All colums except timestamp column
    header = list(raw_file.head(0))[1:]


    # --- Inner loop start
    
    for electrode_index in range(0, len(header)):

        # select electrode index
        #electrode_index = 0

        # Load the selected column signal
        print "Electrode info: Select electrode " + header[electrode_index]
        S = np.array(raw_file[header[electrode_index]])

        # compute length of signal in samples
        L = len(S) # [#samples]
        # compute length of signal in seconds        
        signal_duration = calculate_signal_duration(L, dt)

        # compute power spectral density (PSD) with scipy (DFT using FFT)
        print "Electrode info: Compute PSD"
        f_fft_scipy, Sxx_fft_scipy = signal.periodogram(S, fs=Fs, scaling="density")

        # downsample PSD so that len(P) ~ 5000
        df = f_fft_scipy[1]-f_fft_scipy[0]
        f_downsampled = f_fft_scipy[::int(1/df)]
        P_downsampled = Sxx_fft_scipy[::int(1/df)]

        # select subset of colums using
        # simple low pass filtering using frequency_treshold
        P_downsampled_selected = P_downsampled[0:frequency_treshold+1]
        f_downsampled_selected = f_downsampled[0:frequency_treshold+1]

        # convert elements into string 1.2345e-12 format into 1,2345e-12 because of unscrambler
        P_unscrambler_friendly_list = [ str(e).replace('.',',') for e in P_downsampled_selected ]

        # make pandas series of PDF signal
        P_pandas_s = pandas.Series(P_unscrambler_friendly_list)

        # make pandas dataframe of P_pandas_s with associated electrode number
        P_pandas_df = pandas.DataFrame( {header[electrode_index][0:-5] + " [pV^2/Hz]" : P_pandas_s} ) # removing [pV] part of header string
                                                                                                    # and adding [pV^2/Hz] (Power density in PSD)

        # concatenate the coloumn in P_pandas_df into P_pandas_total_df
        P_pandas_total_df = pandas.concat([P_pandas_total_df, P_pandas_df], axis=1)

        # delete things
        del P_pandas_s
        del P_pandas_df

        # --- inner loop end, go to next electrode_index


    # add other colums with other data and metadata extracted from file_to_import

    # add total energy per frequency column
    E_Hz_list = make_unscrambler_friendly_list([get_total_energy_for_frequency(P_downsampled_selected, int(np.ceil(Hz)), Fs, signal_duration) for Hz in f_downsampled_selected])

    E_Hz_s = pandas.Series(E_Hz_list)
    E_Hz_df = pandas.DataFrame( {"Total energy per Hz [P[Hz]*Fs*duration_seconds for given Hz]" : E_Hz_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, E_Hz_df], axis=1)


    # add signal duration column
    length_time_list = make_unscrambler_friendly_list([signal_duration]*len(P_downsampled_selected))

    l_s = pandas.Series(length_time_list)
    l_df = pandas.DataFrame( {"Duration [s]" : l_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, l_df], axis=1)


    # add total energy of signal column
    E_PSD = calculate_total_signal_energy(P_downsampled, f_downsampled, Fs, signal_duration)
    E_total_list = make_unscrambler_friendly_list([E_PSD]*len(P_downsampled_selected))

    E_total_s = pandas.Series(E_total_list)
    E_total_df = pandas.DataFrame( {"Total energy of signal [sum(P[Hz]*Fs*duration_seconds for all Hz)]" : E_total_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, E_total_df], axis=1)


    # add mean of signal column
    mean_list = make_unscrambler_friendly_list([calculate_signal_mean(S)]*len(P_downsampled_selected))

    mean_s = pandas.Series(mean_list)
    mean_df = pandas.DataFrame( {"Signal mean [pV]" : mean_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, mean_df], axis=1)


    # add min of signal column
    min_list = make_unscrambler_friendly_list([calculate_signal_min(S)]*len(P_downsampled_selected))

    min_s = pandas.Series(min_list)
    min_df = pandas.DataFrame( {"Signal min [pV]" : min_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, min_df], axis=1)


    # add max of signal column
    max_list = make_unscrambler_friendly_list([calculate_signal_max(S)]*len(P_downsampled_selected))

    max_s = pandas.Series(max_list)
    max_df = pandas.DataFrame( {"Signal max [pV]" : max_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, max_df], axis=1)


    # metadata extracted from the string file_to_import

    # add day number of experiment
    day_list = make_unscrambler_friendly_list([days_between(date_seeded, extract_experiment_date(file_to_import))]*len(P_downsampled_selected))

    day_s = pandas.Series(day_list)
    day_df = pandas.DataFrame( {"Day nr. from seeding date " + str(date_seeded) : day_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, day_df], axis=1)


    # add date of experiment yyyy-mm-dd
    date_list = make_unscrambler_friendly_list([extract_experiment_date(file_to_import)]*len(P_downsampled_selected))

    date_s = pandas.Series(date_list)
    date_df = pandas.DataFrame( {"Experiment date [yyyy-mm-dd]" : date_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, date_df], axis=1)


    # add stimulation type
    stimulation_list = make_unscrambler_friendly_list([extract_experiment_stimulation(file_to_import)]*len(P_downsampled_selected))

    stimulation_s = pandas.Series(stimulation_list)
    stimulation_df = pandas.DataFrame( {"Stimulation info" : stimulation_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, stimulation_df], axis=1)

    

    # add conventional frequency bands colums, according to neuroscience
    f_band_list = make_unscrambler_friendly_list([lies_in_frequency_band(int(np.ceil(freq))) for freq in f_downsampled_selected])

    f_band_s = pandas.Series(f_band_list)
    f_band_df = pandas.DataFrame( {"Frequency band" : f_band_s} )

    # nb adding!
    P_pandas_total_df = pandas.concat([P_pandas_total_df, f_band_df], axis=1)


    # lastly, add "Hz" as header of index column
    P_pandas_total_df.index.name = "Hz"
    

    # save P_pandas_total_df to csv
    print "Experiment info: Save to csv"
    P_pandas_total_df.to_csv("PSD/" + file_to_import[0:-4] + "_PSD_pV_and_filename_metadata.csv", sep=";", encoding="ascii")

    # delete things
    del E_Hz_s
    del E_Hz_df
    del l_s
    del l_df
    del E_total_s
    del E_total_df
    del mean_s
    del mean_df
    del min_s
    del min_df
    del max_s
    del max_df
    del day_s
    del day_df
    del date_s
    del date_df
    del stimulation_s
    del stimulation_df
    del f_band_s
    del f_band_df
    
    del P_pandas_total_df
    del raw_file

    # --- outer loop end, go to next experiment_index

# combine all files in PSD folder

# get filenames
filenames = sorted(glob.glob('PSD/*.csv'))

# make empty pandas DataFrame to later save to csv file
P_pandas_combined_df = pandas.DataFrame()

for filename in filenames:
    # Load the csv file 
    print "load " + str(filename)
    raw_file = pandas.read_csv(filename, sep=";", header=0, low_memory=True)

    # concatenate
    print "concatenate"
    P_pandas_combined_df = pandas.concat([P_pandas_combined_df, raw_file], axis=0)
    
    del raw_file

# save P_pandas_total_df to csv
print "save combined csv"
P_pandas_combined_df.to_csv("PSD/combined.csv", sep=";", encoding="ascii", index=False)

del P_pandas_combined_df
