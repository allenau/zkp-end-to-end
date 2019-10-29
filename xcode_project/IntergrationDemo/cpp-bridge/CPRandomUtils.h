// AUTOGENERATED FILE - DO NOT MODIFY!
// This file generated by Djinni from utils.djinni

#import <Foundation/Foundation.h>
@class CPInteger;


@interface CPRandomUtils : NSObject

+ (nonnull NSData *)random:(int32_t)byteLength;

+ (nonnull NSData *)randomWithSeed:(int32_t)byteLength
                              seed:(nonnull NSData *)seed;

+ (nonnull NSString *)randomHex:(int32_t)byteLength;

+ (nonnull NSString *)randomHexWithSeed:(int32_t)byteLength
                                   seed:(nonnull NSData *)seed;

+ (nullable CPInteger *)randomInt:(int32_t)byteLength;

+ (nullable CPInteger *)randomIntWithSeed:(int32_t)byteLength
                                     seed:(nonnull NSData *)seed;

@end
