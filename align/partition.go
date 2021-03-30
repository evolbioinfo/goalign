package align

import (
	"bytes"
	"fmt"
)

type PartitionSet struct {
	names      []string // The name of each partition
	partitions []int    // The partition of each site
	models     []string // The name of the model for each partition
	length     int      // Length of alignment
}

func NewPartitionSet(alignmentLength int) (ps *PartitionSet) {
	partitions := make([]int, alignmentLength)
	for i := 0; i < alignmentLength; i++ {
		partitions[i] = -1
	}
	ps = &PartitionSet{
		names:      make([]string, 0),
		models:     make([]string, 0),
		partitions: partitions,
		length:     alignmentLength,
	}
	return
}

func (ps *PartitionSet) AddRange(partName, modelName string, start, end, modulo int) (err error) {
	if start < 0 {
		err = fmt.Errorf("start of partition is outside of alignment: %d", start)
		return
	}
	if end >= ps.length {
		err = fmt.Errorf("end of partition is outside of alignment: %d", end)
		return
	}
	if modulo <= 0 {
		err = fmt.Errorf("'modulo' value is not authorized: %d", modulo)
		return
	}

	partitionIndex := -1
	for i, p := range ps.names {
		if p == partName {
			partitionIndex = i
			break
		}
	}
	// New partition name
	if partitionIndex == -1 {
		ps.names = append(ps.names, partName)
		ps.models = append(ps.models, modelName)
		partitionIndex = len(ps.names) - 1
	}

	for i := start; i <= end; i += modulo {
		if ps.partitions[i] != -1 {
			err = fmt.Errorf("several partitions are defined for site %d ", i)
			return
		}
		ps.partitions[i] = partitionIndex
	}
	return
}

// If not all sites are in a partition, returns an error
func (ps *PartitionSet) CheckSites() (err error) {
	for j, p := range ps.partitions {
		if p == -1 {
			err = fmt.Errorf("not all sites are in a partition (%d)", j)
			return
		}
	}
	return
}

func (ps *PartitionSet) String() string {
	var buffer bytes.Buffer
	for i, pn := range ps.names {
		buffer.WriteString(ps.models[i])
		buffer.WriteString(",")
		buffer.WriteString(pn)
		buffer.WriteString("=")
		start, end := -1, -1
		first := true
		for j, p := range ps.partitions {

			if p == i {
				if start == -1 {
					start = j
					end = j
				} else if j == end+1 {
					end = j
				}
				if j > end+1 || j == len(ps.partitions)-1 {
					if !first {
						buffer.WriteString(",")
					}
					if start == end {
						buffer.WriteString(fmt.Sprintf("%d", start+1))
					} else {
						buffer.WriteString(fmt.Sprintf("%d-%d", start+1, end+1))
					}
					first = false
					start = j
					end = j
				}
			} else {
				if end != -1 {
					if !first {
						buffer.WriteString(",")
					}
					if start == end {
						buffer.WriteString(fmt.Sprintf("%d", start+1))
					} else {
						buffer.WriteString(fmt.Sprintf("%d-%d", start+1, end+1))
					}
					first = false
				}
				start = -1
				end = -1
			}
		}

		buffer.WriteString("\n")
	}
	return buffer.String()
}

func (ps *PartitionSet) NPartitions() int {
	return len(ps.names)
}

// Returns the partition code associated to the given position
//
// If the position is outside the alignment, then returns -1
func (ps *PartitionSet) Partition(position int) int {
	if position < 0 || position >= len(ps.partitions) {
		return -1
	}

	return ps.partitions[position]
}

// Returns the name of the partition associated to the given index
// If the code does not exist, then returns ""
func (ps *PartitionSet) PartitionName(code int) string {
	if code < 0 || code > len(ps.names) {
		return ""
	}
	return ps.names[code]
}

// Returns the name of the modele associated to the given index
// If the code does not exist, then returns ""
func (ps *PartitionSet) ModeleName(code int) string {
	if code < 0 || code > len(ps.models) {
		return ""
	}
	return ps.models[code]
}

// returns the length of the alignment
func (ps *PartitionSet) AliLength() int {
	return ps.length
}
