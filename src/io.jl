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
	MlemLor

Struct representing a LOR
"""
struct MlemLor
    dx::Float32
    x1::Float32
    y1::Float32
    z1::Float32
    x2::Float32
    y2::Float32
    z2::Float32
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

"""
	write_lors_hdf5(filename, mlor)
	Write lors in a hdf5 format required by petalorust (mlem algo)
"""
function write_lors_hdf5(filename, mlor)

    function set_datatype(::Type{MlemLor})
        dtype = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(MlemLor))
        HDF5.h5t_insert(dtype, "dx", fieldoffset(MlemLor, 1),
                        datatype(Float32))
        HDF5.h5t_insert(dtype, "x1", fieldoffset(MlemLor, 2),
                        datatype(Float32))
        HDF5.h5t_insert(dtype, "y1", fieldoffset(MlemLor, 3),
                        datatype(Float32))
        HDF5.h5t_insert(dtype, "z1", fieldoffset(MlemLor, 4),
                        datatype(Float32))
        HDF5.h5t_insert(dtype, "x2", fieldoffset(MlemLor, 5),
                        datatype(Float32))
        HDF5.h5t_insert(dtype, "y2", fieldoffset(MlemLor, 6),
                        datatype(Float32))
        HDF5.h5t_insert(dtype, "z2", fieldoffset(MlemLor, 7),
                        datatype(Float32))

        HDF5.Datatype(dtype)
    end

    h5f = JPetalo.h5open(filename, "w")

    dtype = set_datatype(MlemLor)
    dspace = dataspace(mlor)
    grp = create_group(h5f, "true_info")
    dset = create_dataset(grp, "lors", dtype, dspace)
    write_dataset(dset, dtype, mlor)


    close(h5f)
end
