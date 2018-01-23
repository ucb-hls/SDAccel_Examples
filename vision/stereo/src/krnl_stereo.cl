/**********
Copyright (c) 2017, Xilinx, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********/

#define FILTER_WIDTH 11
#define FILTER_HEIGHT 11

#define IMAGE_WIDTH 1024
#define IMAGE_HEIGHT 1024

#define B (4)

#define M(x) (((x)-1)/(B) + 1)
#define REG_WIDTH (M(FILTER_WIDTH+B-1)*B)

#if(B == 32)
typedef uint16 bus_t;
#elif(B == 16)
typedef uint8 bus_t;
#elif(B == 8)
typedef uint4 bus_t;
#elif(B == 4)
typedef uint2 bus_t;
#elif(B == 2)
typedef uint bus_t;
#endif

typedef union {
	bus_t b;
	short s[B];
} bus_to_short_t;

void bus_to_short(bus_t in, short out[B]) {
	bus_to_short_t val;

	val.b = in;

	for(int i = 0; i < B; i++) {
		out[i] = val.s[i];
	}
}

bus_t short_to_bus(short in[B]) {
	bus_to_short_t val;

	for(int i = 0; i < B; i++) {
		val.s[i] = in[i];
	}

	return val.b;
}

void get_coef(
	__global bus_t *coef,
	short coef_buf[FILTER_HEIGHT*FILTER_WIDTH]
) {
#ifdef __xilinx__
	__attribute__((xcl_pipeline_loop))
#endif
	for(int i = 0; i < M(FILTER_WIDTH*FILTER_HEIGHT); i++) {
		short tmp[B];
		bus_to_short(coef[i], tmp);

		for(int j = 0; j < B; j++) {
			int k = i*B + j;
			if (k < FILTER_HEIGHT*FILTER_WIDTH) {
				coef_buf[k] = tmp[j];
			}
		}
	}
}

void filter(
	short coef_buf[FILTER_HEIGHT*FILTER_WIDTH],
	global bus_t* input, global bus_t* output
) {

	/* Pad registers to align line_buf read/write */
	short line_reg[FILTER_HEIGHT][REG_WIDTH]
#ifdef __xilinx__
		__attribute__((xcl_array_partition(complete,1)))
		__attribute__((xcl_array_partition(complete,2)))
#endif
		;
	/* Line buffers to store values */
	short line_buf[FILTER_HEIGHT-1][M(IMAGE_WIDTH-REG_WIDTH)*B]
#ifdef __xilinx__
		__attribute__((xcl_array_partition(complete, 1)))
		__attribute__((xcl_array_partition(cyclic, B, 2)))
#endif 
		;

#ifdef __xilinx__
	__attribute__((xcl_pipeline_loop))
#endif
	for(size_t i = 0; i < M(IMAGE_WIDTH*IMAGE_HEIGHT); i++) {
		short input_buf[B];

		/* Read pixels from the input image */
		bus_to_short(input[i], input_buf);

		/* Rotate Buffers */
		for(size_t y = 0; y < FILTER_HEIGHT-1; y++) {
			/* Move the line reg B pixels at a time */
			for(size_t x = 0; x < REG_WIDTH - B; x++) {
				line_reg[y][x] = line_reg[y][x+B];
			}
			/* Add values from line_buf to end of regs */
			for(size_t j = 0; j < B; j++) {
				line_reg[y][(REG_WIDTH - B) + j] = line_buf[y][j + B*(i % (M(IMAGE_WIDTH-REG_WIDTH)))];
			}
			/* Write values from the start of the next line to the line_buf */
			for(size_t j = 0; j < B; j++) {
				line_buf[y][j + B*(i % (M(IMAGE_WIDTH-REG_WIDTH)))] = line_reg[y+1][j];
			}
		}
		/* On last line rotate regs */
		for(size_t x = 0; x < ((M(FILTER_WIDTH+B)-1)*B); x++) {
			line_reg[FILTER_HEIGHT-1][x] = line_reg[FILTER_HEIGHT-1][x+B];
		}
		/* Add the new input data to the end */
		for(size_t j = 0; j < B; j++) {
			line_reg[FILTER_HEIGHT-1][(REG_WIDTH - B) + j] = input_buf[j];
		}

		short filter_sums[B];

		for(size_t j = 0; j < B; j++) {
			int sum = 0;

			for(size_t y = 0; y < FILTER_HEIGHT; y++) {
				for(size_t x = 0; x < FILTER_WIDTH; x++) {
					const size_t offset = REG_WIDTH - FILTER_WIDTH - B + 1;
					short val = line_reg[y][offset + x + j];

					sum += (int) coef_buf[y*FILTER_WIDTH + x] *
					       (int) val;
				}
			}

			/* Handle Saturation */
			if (sum > SHRT_MAX) {
				sum = SHRT_MAX;
			} else if (sum < SHRT_MIN) {
				sum = SHRT_MIN;
			}

			filter_sums[j] = sum;
		}

		output[i] = short_to_bus(filter_sums);
	}
}

//-------------------------------------- TransformPoint ------------------------------------------//
// For efficiency not using union
typedef union {
    float x;
    float y;
}Point2f;

void TransformPoint(const float H[3][3], const float p[2], float ret[2]) {
    float vec_3[3];
    vec_3[0] = p[0];
    vec_3[1] = p[1];
    vec_3[2] = 1;

    float pp[3];
    //matvec(H, vec_3, pp);
    int i, j; 
    for (i = 0; i < 3; i++){
        float tmp;
        for (j = 0; j < 3; j++){
            tmp += H[i][j] * vec_3[j];
        }
        pp[i] = tmp;
    } 
    // float ret[2];
    ret[0] = pp[0]/pp[2];
    ret[1] = pp[1]/pp[2];
}

// JENNY: Matrix3d should be double type 
/*Point2f TransformPoint(const Matrix3d& H, const Point2f& p) {
  // * represents matrix multiplication in opencv
  // # 3 x 3 * 3 x 1
  Vector3d pp = H * Vector3d(p.x, p.y, 1);
  return Point2f(pp.x() / pp.z(), pp.y() / pp.z());
}*/

//-------------------------------------- IsInImage ------------------------------------------//
/*bool IsInImage(const int image_rows, const int image_cols, const float p[2]) {
  if (p[0] >= 0 && p[1] >= 0 &&
      p[0] < image_cols - 0.5 && p[1] < image_rows - 0.5) {
    return true;
  } else {
    return false;
  }
}*/

/*bool IsInImage(const Mat& image, const Point2f& p) {
  if (p.x >= 0 && p.y >= 0 &&
      p.x < image.cols - 0.5 && p.y < image.rows - 0.5) {
    return true;
  } else {
    return false;
  }
}
*/

//-------------------------------------- ExtractPatch ------------------------------------------//
// OpenCV type https://docs.opencv.org/2.4/modules/core/doc/basic_structures.html#mat
// Line 69
#define FirstIfNotNull(X, Y) ( (X) != NULL ? (X) : (Y))

#define WIDTH 100
#define HEIGHT 100
#define NUM_CHAN 1
//template <typename T>
// Assume T is  
// CV_8UC1 unsigned char, with 1 channel 
// CV_8U unsigned  
// Vec3b Vec<unsigned char, 3> 

uint8 ExtractPatch(
                  const unsigned char image [WIDTH][HEIGHT][NUM_CHAN],
                  const float Homography[3][3],
                  const float p[2],
                  const int size[2],
                  int stride,
                  unsigned char patch [WIDTH][HEIGHT][NUM_CHAN], 
                  unsigned char mask [WIDTH][HEIGHT]
                  ) {
  uint8 completed;

  // JENNY see which gets more efficient mapping, not using the opencv library for now 
  //hls::Mat<100, 100, HLS_8UC3> patch();
  //hls::Mat<128, 128, HLS_8UC3> patch(size_width, size_height);
  //hls::Mat<128, 128, HLS_8UC1> size(size_width, size_height);
  // patch->create(size, image.type());
  // *patch = CV_RGB(0, 0, 0);
  // mask->create(size, CV_8U);

  completed = true;
  for (int i = 0; i < WIDTH; ++i) {
    for (int j = 0; j < HEIGHT; ++j) {
      //double p [2];
      float p [2];
      p [0] = (i - (WIDTH + 1) / 2) * stride + p[0]; 
      p [1] = (j - (HEIGHT + 1) / 2) * stride + p[1];
      //double image_point[2];
      float image_point[2];
      TransformPoint(Homography, p, image_point);

      if (p[0] >= 0 && p[1] >= 0 &&
        p[0] < WIDTH - 0.5 && p[1] < HEIGHT - 0.5) {
        for (int k = 0; k < NUM_CHAN; ++k) {
            patch[j][i][k] = image[(int)image_point[0]][(int)image_point[1]][k];
        }
        mask[j][i]= 1;
      } else {
        mask[j][i]= 0;
        completed = false;
      }
    }
  }
}

//}

/*template <typename T>
void ExtractPatch(const cv::Mat& image,
                  const Matrix3d& Homography,
                  const cv::Point2f& p,
                  const cv::Size& size,
                  int stride,
                  Mat *patch_,
                  Mat *mask_ = nullptr,
                  uint8_t *completed_ = nullptr) {
  assert(image.type() == DataType<T>::type);
  Mat local_patch;
  Mat local_mask;
  uint8_t local_completed;
  Mat *patch, *mask;
  uint8_t *completed;
  patch = FirstIfNotNull(patch_, &local_patch);
  mask = FirstIfNotNull(mask_, &local_mask);
  completed = FirstIfNotNull(completed_, &local_completed);
  // Mat patch(size, image.type());
  patch->create(size, image.type());
  *patch = CV_RGB(0, 0, 0);
  mask->create(size, CV_8U);
  *completed = true;
  for (int i = 0; i < size.width; ++i) {
    for (int j = 0; j < size.height; ++j) {
      double x = (i - (size.width + 1) / 2) * stride + p.x;
      double y = (j - (size.height + 1) / 2) * stride + p.y;
      Point2f image_point = TransformPoint(Homography, Point2f(x, y));
      if (IsInImage(image, image_point)) {
        patch->at<T>(j, i) = image.at<T>(image_point);
        mask->at<uint8_t>(j, i) = 1;
      } else {
        mask->at<uint8_t>(j, i) = 0;
        *completed = false;
      }
    }
  }
}

}*/

//-------------------------------------- Homographies ------------------------------------------//

// camera input ?
/*vector<Eigen::Matrix3d> MvsModel::Homographies(double depth) const {
  vector<Matrix3d> homographies;
  for (int i = 0; i < num_cameras(); ++i) {
    homographies.push_back(camera(i)->HomographyFrom(*ref_camera(), depth));
  }
  return homographies;
}*/


//-------------------------------------- CalcNcc ------------------------------------------//
// Only accelerate this function! 
//double CalcNcc(const float p0[WIDTH][HEIGHT], const float p1[WIDTH][HEIGHT]) {
float CalcNcc(const float p0[WIDTH][HEIGHT], const float p1[WIDTH][HEIGHT]) {
  // Why its type is CV_8UC1
  //assert(p0.type() == CV_8UC1);
  float mu0, mu1, sigma0, sigma1, cov;
  mu0 = mu1 = sigma0 = sigma1 = cov = 0;
  int num_pixels = WIDTH * HEIGHT;
  for (int r = 0; r < WIDTH; ++r) {
    for (int c = 0; c < HEIGHT; ++c) {
        mu0 += p0[r][c];
        mu1 += p1[r][c];
    }
  }
  mu0 /= num_pixels;
  mu1 /= num_pixels;

  for (int r = 0; r < WIDTH; ++r) {
    for (int c = 0; c < WIDTH; ++c) {
      sigma0 += (p0[r][c] - mu0) * (p0[r][c] - mu0);
      sigma1 += (p1[r][c] - mu1) * (p1[r][c] - mu1);
      cov += (p0[r][c] - mu0) * (p1[r][c] - mu1);
    }
  }
  return cov / sqrt(sigma0) / sqrt(sigma1);
}

/*double CalcNcc(const Mat& p0, const Mat& p1) {
  assert(p0.cols == p1.cols && p0.rows == p1.rows);
  assert(p0.type() == CV_8UC1);
  double mu0, mu1, sigma0, sigma1, cov;
  mu0 = mu1 = sigma0 = sigma1 = cov = 0;
  int num_pixels = p0.cols * p0.rows;
  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      mu0 += p0.at<uint8_t>(r, c);
      mu1 += p1.at<uint8_t>(r, c);
    }
  }
  mu0 /= num_pixels;
  mu1 /= num_pixels;

  for (int r = 0; r < p0.rows; ++r) {
    for (int c = 0; c < p0.cols; ++c) {
      sigma0 += (p0.at<uint8_t>(r, c) - mu0) * (p0.at<uint8_t>(r, c) - mu0);
      sigma1 += (p1.at<uint8_t>(r, c) - mu1) * (p1.at<uint8_t>(r, c) - mu1);
      cov += (p0.at<uint8_t>(r, c) - mu0) * (p1.at<uint8_t>(r, c) - mu1);
    }
  }
  return cov / sqrt(sigma0) / sqrt(sigma1);
}*/

//-------------------------------------- CalcNccOnDepth ------------------------------------------//
// Line 267

#define NUM_CAMERAS 10
// JENNY: this can only be static 
#define NUM_STEPS 10

int depth_step_ = 10;
int depth_range_[2] = {10, 10};
//float[NUM_CAMERAS][NUM_STEPS][2] CalcNccOnDepth() const {
/*float*** CalcNccOnDepth() {
  float scores[NUM_CAMERAS][NUM_STEPS][2];
  int num_steps = (depth_range_[1] - depth_range_[0]) / depth_step_ + 1;
  // Can make this streaming 
  //for (int i = 0; i < NUM_STEPS; ++i) {
  //  scores[i].resize(num_steps);
  //}
  //ParFor(0, num_steps, [&](int i) {
  int i;
  for(i = 0; i < num_steps; i++ ) {
      float depth = i * depth_step_ + depth_range_[0];
      //JENNY this is not fullly implemented
      float Hs[NUM_CAMERAS][3][3]; 
      //auto Hs = Homographies(depth);
      
      float patches[NUM_STEPS][WIDTH][HEIGHT][3];
      uint8 completed[NUM_CAMERAS];
      for (int j = 0; j < NUM_CAMERAS; ++j) {
        completed[j] = ExtractPatch (
            gray_images_[j], Hs[j],
            ref_point_, patch_size_, patch_stride_,
            &patches[j], nullptr);
      }
      if (!completed[ref_camera_index_]) {
        for (int j = 0; j < NUM_CAMERAS; ++j) {
          scores[j][i][0] = depth;
          scores[j][i][1] = 0.0f;
        }
      } else {
        for (int j = 0; j < NUM_CAMERAS; ++j) {
          if (j == ref_camera_index_) {
            scores[j][i][0] = depth;
            scores[j][i][1] = 1.0f;
          } else if (completed[j]) {
            scores[j][i][0] = depth;
            scores[j][i][1] = CalcNcc(patches[j], patches[ref_camera_index_]));
          } else {
            scores[j][i][0] = depth;
            scores[j][i][1] = 0.0f;
          }
        }
      }
    }
  return scores;
}*/

/*std::vector<std::vector<cv::Point2f>> MvsModel::CalcNccOnDepth() const {
  vector<vector<Point2f>> scores(num_cameras());
  int num_steps = (depth_range_[1] - depth_range_[0]) / depth_step_ + 1;
  // for (int i = 0; i < num_steps; ++i) {
  for (int i = 0; i < num_cameras(); ++i) {
    scores[i].resize(num_steps);
  }
  ParFor(0, num_steps, [&](int i) {
      float depth = i * depth_step_ + depth_range_[0];
      auto Hs = Homographies(depth);
      vector<Mat> patches(num_cameras());
      vector<uint8_t> completed(num_cameras());
      for (int j = 0; j < num_cameras(); ++j) {
        ExtractPatch<uint8_t>(
            gray_images_[j], Hs[j],
            ref_point_, patch_size_, patch_stride_,
            &patches[j], nullptr, &completed[j]);
      }
      if (!completed[ref_camera_index_]) {
        for (int j = 0; j < num_cameras(); ++j) {
          scores[j][i] = Point2f(depth, 0.0f);
        }
      } else {
        for (int j = 0; j < num_cameras(); ++j) {
          if (j == ref_camera_index_) {
            scores[j][i] = Point2f(depth, 1.0f);
          } else if (completed[j]) {
            scores[j][i] = Point2f(
                depth, CalcNcc(patches[j], patches[ref_camera_index_]));
          } else {
            scores[j][i] = Point2f(depth, 0.0f);
          }
        }
      }
    });
  return scores;
}*/

// Line 339
// Runs on the host 
/*void MvsModel::GenerateCostVolume(CostVolume *cost_volume,
                                  DepthMap *depth_map) {
  int num_samples = (depth_range_[1] - depth_range_[0]) / depth_step_;
  cost_volume->Allocate(ref_image_.cols, ref_image_.rows, num_samples);
  (*cost_volume) = 0;
  cost_volume->SetMinDepth(depth_range_[0]);
  cost_volume->SetMaxDepth(depth_range_[1]);
  cost_volume->SetImage(ref_image_);
  depth_map->Allocate(ref_image_.rows, ref_image_.cols);
  int patch_radius = (patch_size_.width - 1) / 2;
  // for (int x = patch_radius; x < ref_image_.cols - patch_radius; ++x) {
  int num_cols = 1000;
  furry::ProgressBar generating_cost_volume("Generating Cost Volume",
                                            // ref_image_.cols - 2 * patch_radius);
                                            num_cols - patch_radius);
  // ParFor(patch_radius, ref_image_.cols - patch_radius, [&](int x) {
  // ParFor(patch_radius, 500, [&](int x) {
  for (int x = patch_radius; x < num_cols; ++x) {
      // LOG(INFO) << "X: " << x;
    for (int y = patch_radius; y < ref_image_.rows - patch_radius; ++y) {
      SetRefPoint(Point2f(x, y));
      auto points = GetPhotoConsistency();
      vector<Point2f> sum_points(points[0].size());
      for (size_t i = 0; i < points[0].size(); ++i) {
        Point2f p(points[0][i].x, 0);
        for (size_t j = 0; j < points.size(); ++j) {
          p.y += points[j][i].y;
        }
        p.y /= points.size();
        sum_points[i] = p;
      }
      float max_label = 0;
      double max_value = sum_points[0].y;
      for (int s = 0; s < num_samples; ++s) {
        cost_volume->Set(y, x, s, 1 - sum_points[s].y);
        if (sum_points[s].y > max_value) {
          max_label = s;
        }
      }
      depth_map->SetDepth(y, x, max_label);
    }
    generating_cost_volume.FinishOne();
    // });
  }
  generating_cost_volume.Done();
}*/


__attribute__((reqd_work_group_size(1,1,1)))
__kernel void krnl_convolve(
   __global bus_t *coef,
   __global bus_t *img_input,
   __global bus_t *img_output
) {
	short coef_buf[FILTER_HEIGHT*FILTER_WIDTH]
#ifdef __xilinx__
		__attribute__((xcl_array_partition(complete, 1)))
#endif
		;

	get_coef(coef, coef_buf);

	filter(coef_buf, img_input, img_output);
}
