#ifndef _IOLIB_CRC32_H_
#define _IOLIB_CRC32_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef IOLIB_CRC
uint32_t iolib_crc32(uint32_t previousCrc32, unsigned char *buf, unsigned int len);
#else
#define iolib_crc32 crc32
#endif

#ifdef __cplusplus
}
#endif

#endif /* _IOLIB_CRC32_H_ */
