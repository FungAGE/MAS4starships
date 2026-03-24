#!/usr/bin/env python3
"""
Setup script for connecting to the original starbase SQLite database.

This script helps you:
1. Set up the SQLite database connection
2. Test the connection
3. Verify that the starbase tables are accessible

Usage:
    python scripts/setup_sqlite_starbase.py
"""

import os
import sys
import django
from django.conf import settings
from django.db import connections
import sqlite3

# Add the project root to Python path
sys.path.append('/mnt/sda/johannesson_lab/adrian/bin/MAS')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'MAS.settings')
django.setup()

def test_sqlite_connection():
    """Test connection to the SQLite database"""
    print("Testing SQLite database connection...")
    
    try:
        # Get the database path from settings
        db_path = settings.DATABASES['starbase']['NAME']
        print(f"Database path: {db_path}")
        
        # Check if file exists
        if not os.path.exists(db_path):
            print(f"❌ Database file not found: {db_path}")
            print("Please ensure the starbase.sqlite file exists at this location")
            return False
        
        # Test connection
        connection = connections['starbase']
        with connection.cursor() as cursor:
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = cursor.fetchall()
            
            print(f"✅ Connected successfully!")
            print(f"Found {len(tables)} tables:")
            for table in tables:
                print(f"  - {table[0]}")
            
            return True
            
    except Exception as e:
        print(f"❌ Connection failed: {e}")
        return False

def inspect_database_schema():
    """Inspect the database schema"""
    print("\nInspecting database schema...")
    
    try:
        connection = connections['starbase']
        with connection.cursor() as cursor:
            # Get all tables
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = [row[0] for row in cursor.fetchall()]
            
            for table in tables:
                print(f"\n--- Table: {table} ---")
                
                # Get table info
                cursor.execute(f"PRAGMA table_info({table})")
                columns = cursor.fetchall()
                
                print("Columns:")
                for col in columns:
                    print(f"  {col[1]} ({col[2]}) - {'NOT NULL' if col[3] else 'NULL'}")
                
                # Get row count
                cursor.execute(f"SELECT COUNT(*) FROM {table}")
                count = cursor.fetchone()[0]
                print(f"Rows: {count}")
                
    except Exception as e:
        print(f"❌ Error inspecting schema: {e}")

def create_sample_models():
    """Create sample model definitions based on the database schema"""
    print("\nCreating sample model definitions...")
    
    try:
        connection = connections['starbase']
        with connection.cursor() as cursor:
            # Get all tables
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = [row[0] for row in cursor.fetchall()]
            
            model_definitions = []
            
            for table in tables:
                # Get table info
                cursor.execute(f"PRAGMA table_info({table})")
                columns = cursor.fetchall()
                
                # Create model class
                model_name = ''.join(word.capitalize() for word in table.split('_'))
                model_def = f"""
class {model_name}(models.Model):
    class Meta:
        managed = False
        db_table = '{table}'
    
"""
                
                # Add fields
                for col in columns:
                    col_name = col[1]
                    col_type = col[2]
                    is_nullable = not col[3]
                    
                    # Map SQLite types to Django field types
                    if col_type.upper() in ['INTEGER', 'INT']:
                        field_type = "models.IntegerField"
                    elif col_type.upper() in ['TEXT', 'VARCHAR']:
                        field_type = "models.TextField"
                    elif col_type.upper() in ['CHAR']:
                        field_type = "models.CharField"
                    else:
                        field_type = "models.TextField"
                    
                    null_param = ", null=True, blank=True" if is_nullable else ""
                    model_def += f"    {col_name} = {field_type}({null_param})\n"
                
                model_def += "\n    def __str__(self):\n        return f\"{model_name} {{self.id}}\"\n"
                model_definitions.append(model_def)
            
            # Write to file
            output_file = os.path.join(settings.BASE_DIR, 'starship', 'starbase_models.py')
            with open(output_file, 'w') as f:
                f.write('"""\n')
                f.write('Starbase models generated from SQLite database schema.\n')
                f.write('These models connect directly to the starbase SQLite database.\n')
                f.write('"""\n\n')
                f.write('from django.db import models\n\n')
                f.write('\n'.join(model_definitions))
            
            print(f"✅ Model definitions written to: {output_file}")
            print("You can now import these models in your views and other parts of the application.")
            
    except Exception as e:
        print(f"❌ Error creating model definitions: {e}")

def main():
    """Main setup function"""
    print("Starbase SQLite Database Setup")
    print("=" * 50)
    
    try:
        # Test connection
        if not test_sqlite_connection():
            print("\n❌ Setup failed: Could not connect to database")
            return
        
        # Inspect schema
        inspect_database_schema()
        
        # Create model definitions
        create_sample_models()
        
        print("\n✅ Setup completed successfully!")
        print("\nNext steps:")
        print("1. Review the generated starbase_models.py file")
        print("2. Update the models as needed for your specific use case")
        print("3. Import the models in your views: from starship.starbase_models import *")
        print("4. Test your application with the new database setup")
        
    except Exception as e:
        print(f"\n❌ Setup failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
