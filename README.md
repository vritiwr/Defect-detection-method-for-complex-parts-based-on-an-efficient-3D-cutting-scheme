# Defect detection on products of complex geometry
This library provides simulation code for An efficient 3D cutting scheme for detecting defects on products of complex geometry.
________________________________________________________________________________________________________
## Requirements
- MATLAB 2018 and newer
## File description
- "slicing.m" is an efficient cutting algorithm.
- "ICP.m" is the iterative closet point  algorithm for the alignment of two 3D models is aligned.
- "defect.m" is used to calculate the defect points after comparing the defective parts with the non-defective parts..
- "Plot_defect_sefment" is used to visualize the resulting defect area..

## Algorithm description
__1. Introduction.__  

```
Surface quality is crucial for high-end equipment functionality, as defects can lead to failures or safety disasters. On-site production inspections are key for ensuring reliability, yet inspecting complex surfaces for hidden defects is costly and often inefficient, with traditional 3D defect detection methods struggling on complex geometries, risking false or missed detections. This study introduces a defect detection scheme using a 3D cutting algorithm for isometric multilayered cutting of parts, aligning each layer's cross-section to detect defects. By converting numerous numerical operations into Boolean operations through reconstructed 3D point cloud data, this approach significantly cuts computing time and mitigates issues from complex geometries. Validated on complex aeroengine impeller parts and other complex structures, our method demonstrates high effectiveness and broad applicability, offering substantial improvements over traditional methods in detecting surface defects.
```

__2. algorithm principle.__  

The pseudo-code for the single-level cutting strategy
```plaintext
// Cutting algorithm and defect detection strategy
𝐈𝐧𝐩𝐮𝐭 𝐃𝐚𝐭𝐚: 1. D = {V, F}
             2. Direction vector (camera system): x(a_1, b_1, c_1), y(a_2, b_2, c_2), z(a_3, b_3, c_3)
             3. Distance parameter \delta z
𝐋𝐨𝐚𝐝𝐢𝐧𝐠 𝐌𝐨𝐝𝐞𝐥s：V, F
𝐒𝐓𝐀𝐆𝐄 𝐈: preprocessing
𝐑𝐨𝐭𝐚𝐭𝐢𝐨𝐧 : Compute 𝑥-axis coordinates
𝐑𝐨𝐭𝐚𝐭𝐢𝐨𝐧 : Compute 𝑦-axis coordinates
𝐓𝐫𝐚𝐧𝐬𝐥𝐚𝐭𝐢𝐨𝐧 : Compute 𝑧-axis coordinates
𝐒𝐓𝐀𝐆𝐄 𝐈𝐈 ∶ Locating intersecting faces
for i = 0:n-1
    if v_iz ≥ 0
       sgn(v_iz) = 1
    else
       sgn(v_iz) = -1
    end if
end for
𝐒𝐓𝐀𝐆𝐄 𝐈𝐈𝐈 ∶ Computing intersection points
for j = 0:m-1
    if|Sgn(f_j)| = 1
      f_new = fj
    end if
    Calculate intersection points by using Eq. (8), Eq. (9)
end for
𝐒𝐓𝐀𝐆𝐄 𝐈𝐕 ∶ Alignment of polygons and defect detection



