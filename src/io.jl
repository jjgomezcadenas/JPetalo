using DataFrames
using HDF5

"""
	PetaloDF

Returns Petalo data sets as DataFrames
"""
struct PetaloDF
	volume_names::Vector{String}
	process_names::Vector{String}
	sensor_xyz::DataFrame
	primaries::DataFrame
	vertices::DataFrame
	total_charge::DataFrame
	waveform::DataFrame
end

"""
	readh5_dset(path::String, folder::String, dset::String)

read an h5 dataset
"""
function readh5_dset(path::String, folder::String, dset::String)
	fid = h5open(path, "r")
	g   = fid[folder]
	dst = read(g, dset)
	return dst
end


"""
	read_abc(path::String)

read abracadabra hdf5 and return relevant data sets
"""
function read_abc(path::String)
	primaries     = DataFrame(readh5_dset(path, "MC", "primaries"))
	sensor_xyz    = DataFrame(readh5_dset(path, "MC", "sensor_xyz"))
	total_charge  = DataFrame(readh5_dset(path, "MC", "total_charge"))
	vertices      = DataFrame(readh5_dset(path, "MC", "vertices"))
	waveform      = DataFrame(readh5_dset(path, "MC", "waveform"))
	volume_names  = readh5_dset(path, "MC", "volume_names")
	process_names = readh5_dset(path, "MC", "process_names")

	primaries[!,"event_id"] = Int.(primaries[!,"event_id"])
	sensor_xyz[!,"sensor_id"] = Int.(sensor_xyz[!,"sensor_id"])
	for col in ["event_id", "sensor_id", "track_id", "parent_id",
		        "process_id", "volume_id", "charge"]

		if col in ["event_id", "sensor_id", "charge"]
			total_charge[!,col] = Int.(total_charge[!,col])
		end

		if col in ["event_id", "sensor_id"]
			waveform[!,col] = Int.(waveform[!,col])
		end

		if col in ["event_id", "track_id", "parent_id","process_id","volume_id"]
			vertices[!,col] = Int.(vertices[!,col])
		end

	end

	return PetaloDF(volume_names,
	                process_names,
					sensor_xyz,
					primaries,
					vertices,
					total_charge,
					waveform)

end
