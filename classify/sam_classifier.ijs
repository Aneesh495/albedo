NB. ==============================================================================
NB. Albedo - Supervised Spectral Angle Mapper (SAM) Classifier
NB.
NB. Classifies every pixel in a hyperspectral data cube into reference classes
NB. by assigning the class with the minimum spectral angle.
NB. Supports maximum angle threshold rejection (label _1 for unclassified/unknown).
NB. ==============================================================================

cocurrent 'albedo_sam_classifier'

load 'detect/sam.ijs'

NB. Classify 3D cube given library matrix M (B x p) and max angle threshold (rad)
NB. Dyadic: (M; max_angle_rad) classify_cube cube
classify_cube=: 4 : 0
  'M max_thresh'=. x
  cube=. y
  'B H W'=. $ cube
  'BM p'=. $ M

  NB. Compute SAM angle map for each class
  angle_maps=. (p, H, W) $ 0.0
  for_c. i. p do.
    target_c=. c {"1 M
    angle_c=. target_c albedo_similarity_sam_cube cube
    angle_maps=. angle_c (<c)} angle_maps
  end.

  NB. For each pixel, find class c with minimum angle
  class_labels=. (H, W) $ _1
  min_angles=. (H, W) $ 0.0

  for_h. i. H do.
    for_w. i. W do.
      pixel_angles=. angle_maps <a:; h; w
      min_ang=. <./ pixel_angles
      min_angles=. min_ang (<h, w)} min_angles
      if. min_ang <: max_thresh do.
        best_c=. (<: +/ pixel_angles = min_ang) { (I. pixel_angles = min_ang)
        class_labels=. best_c (<h, w)} class_labels
      else.
        class_labels=. _1 (<h, w)} class_labels
      end.
    end.
  end.

  class_labels; min_angles
)
