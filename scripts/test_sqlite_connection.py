#!/usr/bin/env python3
"""
Test script to verify SQLite starbase database connection.

This script tests:
1. Database connection
2. Table accessibility
3. Basic queries

Usage:
    python scripts/test_sqlite_connection.py
"""

import os
import sys
import django
from django.conf import settings
from django.db import connections

# Add the project root to Python path
sys.path.append('/mnt/sda/johannesson_lab/adrian/bin/MAS')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'MAS.settings')
django.setup()

def test_database_connection():
    """Test connection to both databases"""
    print("Testing database connections...")
    
    databases = ['default', 'starbase']
    
    for db_name in databases:
        try:
            connection = connections[db_name]
            with connection.cursor() as cursor:
                if db_name == 'starbase':
                    # SQLite specific query
                    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
                    tables = cursor.fetchall()
                    print(f"✅ {db_name} database: Connected successfully")
                    print(f"   Found {len(tables)} tables: {[t[0] for t in tables]}")
                else:
                    # MySQL specific query
                    cursor.execute("SELECT 1")
                    result = cursor.fetchone()
                    print(f"✅ {db_name} database: Connected successfully")
                    
        except Exception as e:
            print(f"❌ {db_name} database: Connection failed - {e}")

def test_starbase_tables():
    """Test access to starbase tables"""
    print("\nTesting starbase table access...")
    
    try:
        connection = connections['starbase']
        with connection.cursor() as cursor:
            # Get all tables
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = [row[0] for row in cursor.fetchall()]
            
            print(f"Found {len(tables)} tables:")
            
            for table in tables:
                try:
                    # Get row count
                    cursor.execute(f"SELECT COUNT(*) FROM {table}")
                    count = cursor.fetchone()[0]
                    print(f"  ✓ {table}: {count} rows")
                    
                    # Get sample data (first row)
                    if count > 0:
                        cursor.execute(f"SELECT * FROM {table} LIMIT 1")
                        sample = cursor.fetchone()
                        print(f"    Sample: {sample[:3]}..." if len(sample) > 3 else f"    Sample: {sample}")
                        
                except Exception as e:
                    print(f"  ✗ {table}: Error - {e}")
                    
    except Exception as e:
        print(f"❌ Error accessing starbase tables: {e}")

def test_cross_database_queries():
    """Test that we can query both databases"""
    print("\nTesting cross-database queries...")
    
    try:
        # Test default database
        default_conn = connections['default']
        with default_conn.cursor() as cursor:
            cursor.execute("SELECT COUNT(*) FROM starship_feature")
            feature_count = cursor.fetchone()[0]
            print(f"✅ Default database: {feature_count} features")
        
        # Test starbase database
        starbase_conn = connections['starbase']
        with starbase_conn.cursor() as cursor:
            # Try to find a starbase table
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'starbase_%' LIMIT 1")
            table = cursor.fetchone()
            
            if table:
                table_name = table[0]
                cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
                count = cursor.fetchone()[0]
                print(f"✅ Starbase database: {count} rows in {table_name}")
            else:
                print("⚠ Starbase database: No starbase_* tables found")
                
    except Exception as e:
        print(f"❌ Error in cross-database queries: {e}")

def main():
    """Main test function"""
    print("SQLite Starbase Database Test")
    print("=" * 50)
    
    try:
        # Test database connections
        test_database_connection()
        
        # Test starbase tables
        test_starbase_tables()
        
        # Test cross-database queries
        test_cross_database_queries()
        
        print("\n" + "=" * 50)
        print("✅ Database connection test completed!")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
