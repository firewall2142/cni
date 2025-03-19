
## Presentation
[![Presentation](https://img.youtube.com/vi/RgtaJ8IN9L0/0.jpg)](https://www.youtube.com/watch?v=RgtaJ8IN9L0)

## Allocation Algorithm

Each day, the information about the number of swabs originating in a district along with location and capacities of labs across the state are provided. An allocation of swabs from collection centers (district collection centers) to testing labs needs to be generated keeping in account the following constraints and cost considerations.

## Allocation Costs

The following are the cost and logistics constraints that need to be considered while generating the allocations:

1.  Cost of testing in a govt. lab is Rs 800 per sample.
2.  Cost of testing in a private lab is Rs 1600 per sample.
3.  Cost of transporting samples outside the source district is equal to Rs 1000 * distance in kilometers from district collection center to centroid of labs. Constraints on the location of labs are given in point 3 in the section below on allocation constraints.
4.  Labs in the same district can be overloaded up to a maximum of 100 samples (over and above a lab's max capacity) for swabs from the same district. Each such overloaded sample will cost Rs 5000 per sample.
5.  Penalty for keeping samples at district headquarters as backlog is Rs 10,000 per sample.

## Allocation Constraints

1.  Sum of allocated swabs to a lab shall not exceed the lab's available capacity for the day (capacity -- backlog), except under the condition given in constraint 2.
2.  A district can overload labs in its same district, up to a maximum of 100 excess swabs per lab. The cost of every such excess swab is given in point 4 of the section on allocation costs.
3.  Allocation of swabs to labs outside of the source district should be such that the maximum distance between any two such labs should be less than 40km.
