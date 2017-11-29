#ifndef INCLUDE_ORA_DATABASECONTAINER_H
#define INCLUDE_ORA_DATABASECONTAINER_H

#include "CondCore/ORA/interface/Handle.h"
#include "RelationalOperation.h"
//
#include <string>
#include <vector>
#include <memory>
#include <typeinfo>

namespace Reflex {
  class Type;
}

namespace ora {

  class IDatabaseSchema;
  class ContainerSchema;
  class DatabaseSession;
  class ContainerUpdateTable;
  class WriteBuffer;
  class UpdateBuffer;
  class ReadBuffer;
  class DeleteBuffer;

  class IteratorBuffer{
    public:
    IteratorBuffer( ContainerSchema& schema, ReadBuffer& buffer );

    ~IteratorBuffer();

    void reset();

    bool next();

    void* getItem();

    void* getItemAsType( const Reflex::Type& type );

    int itemId();

    const Reflex::Type& type();
      
    private:
    SelectOperation m_query;
    int m_itemId;
    ReadBuffer& m_readBuffer;
  };
  
  class DatabaseContainer {
    
    public:
    DatabaseContainer( int contId, const std::string& containerName, const std::string& className,
                       unsigned int containerSize, DatabaseSession& session );

    DatabaseContainer( int contId, const std::string& containerName, const Reflex::Type& containerType,
                       DatabaseSession& session );                  
    
    virtual ~DatabaseContainer();

    int id();

    const std::string& name();

    const std::string& className();

    const Reflex::Type& type();

    const std::string& mappingVersion();

    size_t size();

    void create();

    void drop();

    void extendSchema( const Reflex::Type& dependentType );

    void setAccessPermission( const std::string& principal, bool forWrite );

    Handle<IteratorBuffer> iteratorBuffer();

    bool lock();

    bool isLocked();

    void* fetchItem(int itemId);

    void* fetchItemAsType(int itemId, const Reflex::Type& asType);

    int insertItem( const void* data, const Reflex::Type& type );
    
    void updateItem( int itemId, const void* data, const Reflex::Type& type );    

    void erase( int itemId );
    
    void flush();

    void setItemName( const std::string& name, int itemId );

    bool getNames( std::vector<std::string>& destination );

    private:
    IDatabaseSchema& m_dbSchema;
    std::auto_ptr<ContainerSchema> m_schema;
    std::auto_ptr<WriteBuffer> m_writeBuffer;
    std::auto_ptr<UpdateBuffer> m_updateBuffer;
    std::auto_ptr<ReadBuffer>   m_readBuffer;
    std::auto_ptr<DeleteBuffer> m_deleteBuffer;
    Handle<IteratorBuffer> m_iteratorBuffer;
    size_t m_size;
    ContainerUpdateTable& m_containerUpdateTable;
    bool m_lock;
  };

}

#endif
