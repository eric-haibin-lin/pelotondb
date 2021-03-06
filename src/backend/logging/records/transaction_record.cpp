
#include <iostream>

#include "backend/logging/records/transaction_record.h"

namespace peloton {
namespace logging {

/**
 * @brief Serialize given data
 * @return true if we serialize data otherwise false
 */
bool TransactionRecord::Serialize(CopySerializeOutput &output) {
  bool status = true;
  output.Reset();

  // First, write out the log record type
  output.WriteEnumInSingleByte(log_record_type);

  // Then reserve 4 bytes for the header size to be written later
  size_t start = output.Position();
  output.WriteInt(0);
  output.WriteLong(txn_id);

  // Write out the header now
  int32_t header_length =
      static_cast<int32_t>(output.Position() - start - sizeof(int32_t));
  output.WriteIntAt(start, header_length);

  message_length = output.Size();
  message = new char[message_length];
  memcpy(message, output.Data(), message_length);

  return status;
}

/**
 * @brief Deserialize LogRecordHeader
 * @param input
 */
void TransactionRecord::Deserialize(CopySerializeInputBE &input) {
  // Get the message length
  input.ReadInt();

  // Just grab the transaction id
  txn_id = (txn_id_t)(input.ReadLong());
}

// Used for peloton logging
size_t TransactionRecord::GetTransactionRecordSize(void) {
  // log_record_type + header_legnth + transaction_id
  return sizeof(char) + sizeof(int) + sizeof(long);
}

void TransactionRecord::Print(void) {
  std::cout << "#LOG TYPE:" << LogRecordTypeToString(GetType()) << "\n";
  std::cout << " #Txn ID:" << GetTransactionId() << "\n";
  std::cout << "\n";
}

}  // namespace logging
}  // namespace peloton
