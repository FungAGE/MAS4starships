class DatabaseRouter:    
    """
    Database router to direct starbase models to the SQLite starbase database
    and all other models to the default MySQL database.
    """
    
    # Models that should use the starbase database (SQLite)
    starbase_models = {
        'accessions', 'ships', 'captains', 'taxonomy', 'genome', 'papers',
        'familynames', 'starshipfeatures', 'navis', 'haplotype', 'gff', 'joinedships'
    }
    
    def db_for_read(self, model, **hints):
        # Check if model is in starbase_models or has starbase_ prefix in table name
        if (model._meta.model_name in self.starbase_models or 
            model._meta.db_table.startswith('starbase_')):
            return 'starbase'
        return 'default'

    def db_for_write(self, model, **hints):
        # Check if model is in starbase_models or has starbase_ prefix in table name
        if (model._meta.model_name in self.starbase_models or 
            model._meta.db_table.startswith('starbase_')):
            return 'starbase'
        return 'default'

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations between models in the same database.
        For cross-database relations, we need to be more careful.
        """
        db_set = {obj1._state.db, obj2._state.db}
        if len(db_set) == 1:
            # Same database
            return True
        else:
            # Cross-database relations - allow but be cautious
            # Note: Some Django features may not work across databases
            return True

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        Ensure starbase models only migrate to starbase database
        and other models only migrate to default database.
        """
        if db == 'starbase':
            # Only allow starbase models to migrate to starbase database
            return model_name in self.starbase_models
        elif db == 'default':
            # Only allow non-starbase models to migrate to default database
            return model_name not in self.starbase_models
        return False
