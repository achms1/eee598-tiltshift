#include <jni.h>
#include <arm_neon.h>
#include <iostream>
#include <string>
#include <cpu-features.h>
#include <cmath>

//for logging
#include <android/log.h>

#define LOG_TAG "testjni"
#define PI 3.14

using namespace std;
//__android_log_print(ANDROID_LOG_INFO, "TAG", "%d, %d", width, height);
// Functions declared for calculating Gaussian Kernal, Intermediate and Final matrix
float * create_gaussian_kernel(int, float);
void copy_row(jint *, jint *, int, jint);
int pvector(jint *, int, int, int, float *, int, int);
int qvector(jint *, int, int, int, float *, int, int);

extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftcppnative(JNIEnv *env,
                                                                               jclass instance,
                                                                               jintArray inputPixels_,
                                                                               jintArray outputPixels_,
                                                                               jint width,
                                                                               jint height,
                                                                               jfloat sigma_far,
                                                                               jfloat sigma_near,
                                                                               jint a0, jint a1,
                                                                               jint a2, jint a3) {
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "sigma far = %f, sigma near = %f", sigma_far, sigma_near);

    long length = env->GetArrayLength(inputPixels_);
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "width = %d, length = %d", width, height);
    float *kern = NULL;
    float sigma=0.0;
    int r_prev = 0, r_new = 0;

    for (int j=0;j<height;j++){

        //calculate sigma value based on the pixel position
        if(j<a0){
            sigma = sigma_far;
        }
        else if(j<a1){
            sigma = (sigma_far*(a1-j))/(a1-a0);
        }
        else if(j>=a1 && j<a2){
            sigma = 0.5;
        }
        else if(j<a3){
            sigma = (sigma_near*(j-a2))/(a3-a2);
        }
        else if(j<height){
            sigma = sigma_near;
        }

        //__android_log_print(ANDROID_LOG_INFO, "TAG", "sigma = %f", sigma);
        if(sigma<0.6){
            // Calling in copy_row function to copy the pixel values when sigma < 0.6 (Gaussian blur should not be applied)
            copy_row(pixels, outputPixels, j, width);
            continue;
        }

        r_new= int(ceil(2*sigma));
        //__android_log_print(ANDROID_LOG_INFO, "TAG", "height = %d, current = %d", height, j);

        if(r_new!=r_prev){
            delete[] kern;
            kern = create_gaussian_kernel(r_new, sigma);
            r_prev = r_new;
        }

        for (int i=0;i<width;i++) {


            int p = pvector(pixels, i, j, r_new, kern, height, width);

            //int color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 8 | (B & 0xff);
            outputPixels[j*width+i]= p;
        }
    }
    delete kern;
    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    __android_log_print(ANDROID_LOG_INFO, "TAG", "operation finished!\n");
    return 0;
}

// Function to calculate intermediate vector using weight vector approach(fast approach)
int qvector(int * pixels, int x, int y, int r, float * kern, int height, int width) {
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "inside q");
    int len = (2*r)+1;
    int q;

    int tempB = 0, tempR = 0, tempG = 0;
    for(int i=0; i<len; i++){
        if( (y-r+i) < 0 || (y-r+i) >= height){
            q = 0;
        }else {
            q = pixels[(y - r + i) * width + x];
        }
        int B = q&0xff;
        int G = (q>>8)&0xff;
        int R = (q>>16)&0xff;

        tempB += B*kern[i];
        tempG += G*kern[i];
        tempR += R*kern[i];
    }
    int A = 0xff;
    int colorq = (A & 0xff) << 24 | (tempR & 0xff) << 16 | (tempG & 0xff) << 8 | (tempB & 0xff);
    return colorq;
}
// Function to calculate the final output vector using weight vector approach(fast approach)
int pvector(int * pixels, int x, int y, int r, float * kern, int height, int width){
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "inside p1");
    int len = (2*r)+1;
    int p,q;

    float tempB = 0, tempR = 0, tempG = 0;
    for(int i=0; i<len; i++){
        if( (x-r+i) < 0 || (x-r+i) >= width){
            //__android_log_print(ANDROID_LOG_INFO, "TAG", "inside p2");
            q = 0;
        }else {
            // Calling in qvector function to calculate intermediate vector
            q = qvector(pixels, x-r+i, y, r, kern, height, width);
        }
        //__android_log_print(ANDROID_LOG_INFO, "TAG", "inside p3");
        int B = q&0xff;
        int G = (q>>8)&0xff;
        int R = (q>>16)&0xff;
        //__android_log_print(ANDROID_LOG_INFO, "TAG", "inside p4");
        tempB += (B*kern[i]);
        tempG += (G*kern[i]);
        tempR += (R*kern[i]);
        //__android_log_print(ANDROID_LOG_INFO, "TAG", "inside p5");
    }
    int A = 0xff;
    p = (A & 0xff) << 24 | (int(tempR) & 0xff) << 16 | (int(tempG) & 0xff) << 8 | (int(tempB) & 0xff);
    return p;
}

// Function to copy the pixel values when sigma < 0.6
void copy_row(jint * pixels, jint * output, int j, jint width){
    for(int i=0; i<width; i++){
        int B = pixels[j*width+i]&0xff;
        int G = (pixels[j*width+i]>>8)&0xff;
        int R = (pixels[j*width+i]>>16)&0xff;
        int A = 0xff;

        int color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 8 | (B & 0xff);
        output[j*width+i]=color;
    }
}


// Function to calculate Gaussian Kernal matrix
float * create_gaussian_kernel(int r, float sigma){
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "function called!");
    int len = (2*r)+1;
    float * G = new float [len];
    for(int i=0; i<len; i++){
        float temp = -1*((-1*r)+i)*((-1*r)+i);
        temp /= (2*sigma*sigma);
        temp = exp(temp);
        temp /= sqrt(2*PI*sigma*sigma);
        G[i] = (temp);
        //__android_log_print(ANDROID_LOG_INFO, "TAG", "kern[%d] = %f", i, G[i]);
    }
    return G;
}

void print_func(uint16x8_t data){
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "yo");
    uint16_t p[8];

    vst1q_u16(p, data);

    for(int i=0; i<8; i++){
        __android_log_print(ANDROID_LOG_INFO, "TAG", "element[%d] : %d\n", i, p[i]);
    }
}

// Function to calculate Gaussian Kernal matrix
int * neon_create_gaussian_kernel(int r, float sigma){
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "function called!");
    int len = (2*r)+1;
    int * G = new int [len];
    for(int i=0; i<len; i++){
        float temp = -1*((-1*r)+i)*((-1*r)+i);
        temp /= (2*sigma*sigma);
        temp = exp(temp);
        temp /= sqrt(2*PI*sigma*sigma);
        temp *= 128;
        G[i] = (int)(ceil(temp));
        //__android_log_print(ANDROID_LOG_INFO, "TAG", "kern[%d] = %d", i, G[i]);
    }
    return G;
}

static int16x8x4_t temp_sum; //accumulator
static uint16x8_t temp_product;
static uint16x8_t big_blue;
static uint16x8_t big_green;
static uint16x8_t big_red;

static uint8x8x4_t pix;
static uint8x8x4_t tempgk;
static uint8x8_t gk;
static uint16x8_t big_gk;

int neon_dotproduct(int * vector1, int * vector2, int len) {
    int vectorSize = 8*4;
    int parts = len / 8;
    int remainder = len % 8;

    //__android_log_print(ANDROID_LOG_INFO, "TAG", "len = %d, parts = %d, remainder = %d", len, parts, remainder);
    temp_sum.val[0] = vdupq_n_u16(0);
    temp_sum.val[1] = vdupq_n_u16(0);
    temp_sum.val[2] = vdupq_n_u16(0);
    temp_sum.val[3] = vdupq_n_u16(0);

    uint8_t * pA = (uint8_t *)vector1;
    uint8_t *pB = (uint8_t *)vector2;

    for(int i = 0; i < parts; i++) {
        pix = vld4_u8(pA + (i*vectorSize));
        tempgk = vld4_u8(pB + (i*vectorSize));
        gk = tempgk.val[0];
        big_gk = vmovl_u8(gk);

        big_blue = vmovl_u8(pix.val[0]);
        big_green = vmovl_u8(pix.val[1]);
        big_red = vmovl_u8(pix.val[2]);

        temp_product = vmulq_u16(big_blue, big_gk);
        temp_product = vshrq_n_u16(temp_product, 7);
        temp_sum.val[0] = vaddq_u16(temp_sum.val[0], temp_product);

        temp_product = vmulq_u16(big_green, big_gk);
        temp_product = vshrq_n_u16(temp_product, 7);
        temp_sum.val[1] = vaddq_u16(temp_sum.val[1], temp_product);

        temp_product = vmulq_u16(big_red, big_gk);
        temp_product = vshrq_n_u16(temp_product, 7);
        temp_sum.val[2] = vaddq_u16(temp_sum.val[2], temp_product);

    }

    if(remainder!=0){
        //__android_log_print(ANDROID_LOG_INFO, "TAG", "Remainder\n");
        pix = vld4_u8(pA + (parts*vectorSize));
        tempgk = vld4_u8(pB + (parts*vectorSize));
        gk = tempgk.val[0];
        big_gk = vmovl_u8(gk);

        big_blue = vmovl_u8(pix.val[0]);
        big_green = vmovl_u8(pix.val[1]);
        big_red = vmovl_u8(pix.val[2]);

        temp_product = vmulq_u16(big_blue, big_gk);
        temp_product = vshrq_n_u16(temp_product, 7);
        temp_sum.val[0] = vaddq_u16(temp_sum.val[0], temp_product);

        temp_product = vmulq_u16(big_green, big_gk);
        temp_product = vshrq_n_u16(temp_product, 7);
        temp_sum.val[1] = vaddq_u16(temp_sum.val[1], temp_product);

        temp_product = vmulq_u16(big_red, big_gk);
        temp_product = vshrq_n_u16(temp_product, 7);
        temp_sum.val[2] = vaddq_u16(temp_sum.val[2], temp_product);
    }

    int result=0xff;
    static unsigned short channel[8];

    unsigned int channel_sum = 0;
    vst1q_u16(channel, temp_sum.val[2]);
    for (int l = 0; l < 8; l++) {
        channel_sum += channel[l];
    }
    if(channel_sum > 255)
        channel_sum = 255;
    result = result<<8;
    result += channel_sum;

    channel_sum = 0;
    vst1q_u16(channel, temp_sum.val[1]);
    for (int l = 0; l < 8; l++) {
        channel_sum += channel[l];
    }
    if(channel_sum > 255)
        channel_sum = 255;
    result = result<<8;
    result += channel_sum;

    channel_sum = 0;
    vst1q_u16(channel, temp_sum.val[0]);
    for (int l = 0; l < 8; l++) {
        channel_sum += channel[l];
    }
    if(channel_sum > 255)
        channel_sum = 255;
    result = result<<8;
    result += channel_sum;

    return result;

}

int neon_qvector(int * pixels, int x, int y, int r, int * kern, int height, int width) {
    int len = (2*r)+1;

    int * p = NULL;
    p = new int [len];
    for(int i=0; i<len; i++){
        if( (y-r+i) < 0 || (y-r+i) >= height){
            p[i] = 0;
        }
        else {
            p[i] = pixels[(y - r + i) * width + x];
        }
    }

    int res = neon_dotproduct(p, kern, len);
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "result=%d\n", res);
    delete p;
    return res;
}

void neon_pvector(int * pixels, int * outputPixels, int x, int y, int r, int * kern, int height, int width){
    int len = (2*r)+1;
    int * qvector = new int [len];

    for(int i=0; i<len; i++){
        if( (x-r+i) < 0 || (x-r+i) >= width){
            qvector[i] = 0;
        }
        else {
            qvector[i] = neon_qvector(pixels, x - r + i, y, r, kern, height, width);
        }
    }

    int result = neon_dotproduct(qvector, kern, len);
    //__android_log_print(ANDROID_LOG_INFO, "TAG", "pixels = %d\n", result);
    outputPixels[y*width + x] = result;
    delete[] qvector;

    return;
}

extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftneonnative(JNIEnv *env,
                                                                                jclass instance,
                                                                                jintArray inputPixels_,
                                                                                jintArray outputPixels_,
                                                                                jint width,
                                                                                jint height,
                                                                                jfloat sigma_far,
                                                                                jfloat sigma_near,
                                                                                jint a0, jint a1,
                                                                                jint a2, jint a3) {
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);

    //__android_log_print(ANDROID_LOG_INFO, "TAG", "width = %d, length = %d", width, height);
    int *kern = NULL;
    float sigma=0.0;
    int r_prev = 0, r_new = 0;

    for (int j=0;j<height;j++){
        if(j<a0){
            sigma = sigma_far;
        }
        else if(j<a1){
            sigma = (sigma_far*(a1-j))/(a1-a0);
        }
        else if(j>=a1 && j<a2){
            sigma = 0.5;
        }
        else if(j<a3){
            sigma = (sigma_near*(j-a2))/(a3-a2);
        }
        else if(j<height){
            sigma = sigma_near;
        }

        //__android_log_print(ANDROID_LOG_INFO, "TAG", "sigma = %f", sigma);
        if(sigma<0.6){
            // Calling in copy_row function to copy the pixel values when sigma < 0.6 (Gaussian blur should not be applied)
            copy_row(pixels, outputPixels, j, width);
            continue;
        }

        r_new= int(ceil(2*sigma));
        //__android_log_print(ANDROID_LOG_INFO, "TAG", "height = %d, current = %d", height, j);

        if(r_new!=r_prev){
            delete[] kern;
            kern = neon_create_gaussian_kernel(r_new, sigma);
            r_prev = r_new;
        }

        for (int i=0;i<width;i++) {
           //if(i==10&&j==115)
                neon_pvector(pixels, outputPixels, i, j, r_prev, kern, height, width);
        }
    }

    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    __android_log_print(ANDROID_LOG_INFO, "TAG", "operation finished!\n");
    return 0;

}