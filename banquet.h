/** Parameter set names */
typedef enum {
  PARAMETER_SET_INVALID = 0,
  Banquet_L1_Param1 = 1,
  PARAMETER_SET_MAX_INDEX = 2
} banquet_params_t;

/** Public key */
typedef struct {
  uint8_t data[PICNIC_MAX_PUBLICKEY_SIZE];
} picnic_publickey_t;

/** Private key */
typedef struct {
  uint8_t data[PICNIC_MAX_PRIVATEKEY_SIZE];
} picnic_privatekey_t;